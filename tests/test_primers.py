"""Tests for primer design."""
import pytest
from gsynth_core.primers import (
    design_primer_pair, design_cloning_primers, analyze_primer, PrimerParams,
)


TEMPLATE = (
    "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCA"
    "GCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAA"
    "TGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTC"
)


class TestPrimerDesign:
    def test_basic_pair(self):
        pair = design_primer_pair(TEMPLATE)
        assert pair is not None
        assert pair.tm_difference < 3.0
        # Forward primer should match the start of the template
        assert TEMPLATE.upper().find(pair.forward.sequence) >= 0

    def test_cloning_primer_sites(self):
        pair = design_cloning_primers(TEMPLATE, left_enzyme="EcoRI", right_enzyme="BamHI")
        assert pair is not None
        assert "GAATTC" in pair.forward.sequence
        assert "GGATCC" in pair.reverse.sequence

    def test_no_solution_strict_params(self):
        # Force impossible constraint
        params = PrimerParams(min_tm=95.0, max_tm=99.0)
        pair = design_primer_pair(TEMPLATE, params=params)
        assert pair is None


class TestPrimerAnalysis:
    def test_analyze_returns_full_report(self):
        an = analyze_primer("CCATGATTACGGATTCACTGG")
        assert an.length == 21
        assert an.tm > 0
        assert 0 <= an.gc_content <= 100
        assert an.hairpin.risk_level in ("low", "medium", "high")
