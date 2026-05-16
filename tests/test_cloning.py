"""Tests for cloning modules."""
import pytest
from gsynth_core.cloning import (
    design_gibson_overlaps, simulate_gibson_assembly,
    design_golden_gate, simulate_golden_gate_assembly,
    check_ligation_compatibility, simulate_restriction_cloning,
)
from gsynth_core.restriction.enzymes import OverhangType


class TestGibson:
    def test_overlap_design(self):
        frags = [("a", "A" * 50 + "GCATGCATGCATGCA"),
                 ("b", "GCATGCATGCATGCA" + "T" * 50)]
        designed = design_gibson_overlaps(frags, circular=False, min_overlap_len=10)
        assert designed[0].right_overlap
        # Overlap should equal the shared region
        assert designed[0].right_overlap == designed[1].left_overlap

    def test_min_fragments_validation(self):
        with pytest.raises(ValueError):
            design_gibson_overlaps([("only", "ATGC")])


class TestGoldenGate:
    def test_assembly_removes_enzyme_site(self):
        parts = [("p1", "ATGAAACGCCTGGCTGTT"), ("p2", "ATGGTGAGCAAGGGC"),
                 ("p3", "GACTACAAGGACGACG"), ("p4", "TAATAAGCGCGCGCG")]
        gg = design_golden_gate(parts, enzyme="BsaI")
        assembled = simulate_golden_gate_assembly(gg, circular=True)
        # BsaI site (GGTCTC) and its RC (GAGACC) must not appear
        assert "GGTCTC" not in assembled
        assert "GAGACC" not in assembled

    def test_assembly_preserves_inserts(self):
        parts = [("p1", "AAAAAAAAAAA"), ("p2", "TTTTTTTTTTT"),
                 ("p3", "GGGGGGGGGGG"), ("p4", "CCCCCCCCCCC")]
        gg = design_golden_gate(parts, enzyme="BsaI")
        assembled = simulate_golden_gate_assembly(gg, circular=True)
        # All four sequences must appear (not necessarily their flanking 4-nt overhangs)
        for _, ins in parts:
            assert ins in assembled


class TestLigation:
    def test_compatible_overhangs(self):
        assert check_ligation_compatibility(
            "AATT", OverhangType.FIVE_PRIME, "AATT", OverhangType.FIVE_PRIME,
        )

    def test_incompatible_overhangs(self):
        assert not check_ligation_compatibility(
            "AATT", OverhangType.FIVE_PRIME, "GATC", OverhangType.FIVE_PRIME,
        )

    def test_blunt_blunt(self):
        assert check_ligation_compatibility(
            "", OverhangType.BLUNT, "", OverhangType.BLUNT,
        )

    def test_blunt_and_sticky_dont_mix(self):
        assert not check_ligation_compatibility(
            "AATT", OverhangType.FIVE_PRIME, "", OverhangType.BLUNT,
        )
