"""Tests for restriction enzymes and digestion."""
import pytest
from gsynth_core.restriction import (
    find_sites, get_enzyme, list_enzymes, simulate_digestion,
    simulate_double_digestion, suggest_compatible_ends, OverhangType,
)


class TestEnzymes:
    def test_get_enzyme(self):
        e = get_enzyme("EcoRI")
        assert e.site == "GAATTC"
        assert e.overhang_sequence == "AATT"
        assert e.overhang_type == OverhangType.FIVE_PRIME

    def test_blunt_cutter(self):
        e = get_enzyme("SmaI")
        assert e.overhang_type == OverhangType.BLUNT

    def test_three_prime_overhang(self):
        e = get_enzyme("PstI")
        assert e.overhang_type == OverhangType.THREE_PRIME

    def test_palindromic(self):
        assert get_enzyme("EcoRI").is_palindromic
        assert get_enzyme("BamHI").is_palindromic

    def test_list_enzymes(self):
        names = list_enzymes(source="curated")
        assert len(names) > 70
        assert "EcoRI" in names

    def test_compatible_ends(self):
        # BamHI/BglII/BclI all have GATC overhang
        compat = suggest_compatible_ends("BamHI")
        assert "BglII" in compat
        assert "BclI" in compat

    def test_unknown(self):
        with pytest.raises(KeyError):
            get_enzyme("FakeI_999")


class TestFindSites:
    def test_finds_ecori(self):
        hits = find_sites("AAAAGAATTCAAAA")
        names = [h.enzyme.name for h in hits]
        assert "EcoRI" in names

    def test_iupac_site(self):
        # ApoI = RAATTY → should hit AAATTC and GAATTC
        hits = find_sites("AAATTCAAAGAATTCAAA", enzymes=["ApoI"])
        assert len(hits) >= 2


class TestDigest:
    def test_single_cut(self):
        frags = simulate_digestion("AAAAGAATTCAAAA", "EcoRI")
        assert len(frags) == 2

    def test_double_digest(self):
        frags = simulate_double_digestion(
            "AAAAGAATTCAAAAGGATCCAAAAA", "EcoRI", "BamHI",
        )
        assert len(frags) == 3
        # check overhangs
        assert frags[0].right_overhang == "AATT"
        assert frags[1].left_overhang == "AATT"
        assert frags[1].right_overhang == "GATC"

    def test_no_cut(self):
        frags = simulate_digestion("AAAAAAA", "EcoRI")
        assert len(frags) == 1
        assert frags[0].length == 7

    def test_circular(self):
        # plasmid with 2 EcoRI sites
        seq = "GAATTC" + "AAAAAAAAAA" + "GAATTC" + "AAAAAAAAAA"
        frags = simulate_digestion(seq, "EcoRI", circular=True)
        assert len(frags) == 2
