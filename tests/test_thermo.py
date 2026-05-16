"""Tests for Tm calculations."""
import pytest
from gsynth_core.thermo import (
    melting_temperature, melting_temperature_nn,
    melting_temperature_marmur, melting_temperature_gc, TmMethod,
)


class TestTm:
    def test_m13_universal(self):
        """M13 universal primer ~56°C (literature)."""
        tm = melting_temperature("GTAAAACGACGGCCAGT")
        assert 54 < tm < 62, f"M13 Tm out of range: {tm}"

    def test_short_primer(self):
        tm = melting_temperature("ATGCATGC")
        assert 0 < tm < 40

    def test_long_primer(self):
        tm = melting_temperature("ATGCATGCATGCATGCATGCATGCATGC")
        assert 50 < tm < 80

    def test_marmur_short(self):
        tm = melting_temperature_marmur("ATGC")
        # 2*2 + 4*2 = 12
        assert tm == 12.0

    def test_gc_method(self):
        tm = melting_temperature_gc("ATGCATGCATGCATGC")
        assert tm > 0

    def test_invalid_input(self):
        with pytest.raises(ValueError):
            melting_temperature("XYZ")

    def test_method_marmur(self):
        tm = melting_temperature("ATGCATGCATGC", method=TmMethod.MARMUR)
        # A+T = 6, G+C = 6 → 2*6 + 4*6 = 36 (no penalty for len < 13)
        assert tm == 36.0

    def test_palindrome_lower_tm(self):
        """Palindromic sequences should have lower Tm than non-palindromic."""
        # Both have same length and GC%
        pal = melting_temperature_nn("GAATTCGAATTC")
        non_pal = melting_temperature_nn("GAATTCAAGCTT")  # not palindromic
        assert pal < non_pal  # symmetry penalty
