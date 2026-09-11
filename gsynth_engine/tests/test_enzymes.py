import pytest

from gsynth_engine.constants import (
    ALL_ENZYMES,
    RESTRICTION_ENZYMES,
    left_remainders,
    overhang,
    right_remainders,
    supplies_start_codon,
)
from gsynth_engine.sequence import reverse_complement

biopython = pytest.importorskip("Bio.Restriction", reason="regeneration source")


class TestTheCuratedSetIsCorrect:


    @pytest.mark.parametrize("name", sorted(RESTRICTION_ENZYMES))
    def test_geometry_matches_rebase(self, name):

        enzyme = getattr(biopython, name)
        site = str(enzyme.site)
        ours = RESTRICTION_ENZYMES[name]

        assert ours["recognition"] == site
        assert ours["cut_top"] == enzyme.fst5
        assert ours["cut_bottom"] == len(site) + enzyme.fst3


class TestTheWideTable:
    def test_contains_every_curated_enzyme(self):

        missing = sorted(set(RESTRICTION_ENZYMES) - set(ALL_ENZYMES))
        assert missing == []
        for name, spec in RESTRICTION_ENZYMES.items():
            for key in ("recognition", "cut_top", "cut_bottom"):
                assert ALL_ENZYMES[name][key] == spec[key], name

    def test_is_substantially_wider_than_the_freezer(self):

        assert len(ALL_ENZYMES) >= 100

    @pytest.mark.parametrize("name", sorted(ALL_ENZYMES))
    def test_every_entry_is_a_real_type_IIP_specification(self, name):

        spec = ALL_ENZYMES[name]
        site = str(spec["recognition"])
        top, bottom = int(spec["cut_top"]), int(spec["cut_bottom"])

        assert site and not set(site) - set("ACGT"), site
        assert 0 <= top <= len(site)
        assert 0 <= bottom <= len(site)

    @pytest.mark.parametrize("name", sorted(ALL_ENZYMES))
    def test_the_derived_ends_reassemble_the_site(self, name):

        site = str(ALL_ENZYMES[name]["recognition"])
        left_fwd, left_rev = left_remainders(name)
        right_fwd, right_rev = right_remainders(name)

        assert right_fwd + left_fwd == site
        assert left_rev + right_rev == reverse_complement(site)

    @pytest.mark.parametrize("name", sorted(ALL_ENZYMES))
    def test_overhang_agrees_with_the_cut_positions(self, name):

        spec = ALL_ENZYMES[name]
        site = str(spec["recognition"])
        top, bottom = int(spec["cut_top"]), int(spec["cut_bottom"])
        sequence, kind = overhang(name)

        if top == bottom:
            assert (sequence, kind) == ("", "blunt")
        else:
            low, high = min(top, bottom), max(top, bottom)
            assert sequence == site[low:high]
            assert kind == ("5'" if top < bottom else "3'")

    def test_start_codon_supply_depends_on_the_retained_remainder(self):
        assert {name for name in ALL_ENZYMES if supplies_start_codon(name)} == {
            "CviAII", "FatI", "NdeI",
        }


        for name in ("CciI", "FaeI", "NcoI", "NsiI", "PciI", "SphI"):
            assert "ATG" in str(ALL_ENZYMES[name]["recognition"])
            assert not supplies_start_codon(name)


class TestIsoschizomersAreCollapsed:
    def test_one_name_per_cut_specification(self):

        seen: dict[tuple, str] = {}
        for name, spec in ALL_ENZYMES.items():
            key = (spec["recognition"], spec["cut_top"], spec["cut_bottom"])
            assert key not in seen, f"{name} duplicates {seen.get(key)}"
            seen[key] = name

    def test_the_familiar_name_survives(self):

        assert ALL_ENZYMES["NheI"]["recognition"] == "GCTAGC"
        assert "AsuNHI" not in ALL_ENZYMES, "a true isoschizomer of NheI"

    def test_a_neoschizomer_is_not_collapsed_into_its_twin(self):

        nhe, bmt = ALL_ENZYMES["NheI"], ALL_ENZYMES["BmtI"]
        assert nhe["recognition"] == bmt["recognition"] == "GCTAGC"
        assert (nhe["cut_top"], nhe["cut_bottom"]) == (1, 5)
        assert (bmt["cut_top"], bmt["cut_bottom"]) == (5, 1)
        assert overhang("NheI") == ("CTAG", "5'")
        assert overhang("BmtI") == ("CTAG", "3'")
