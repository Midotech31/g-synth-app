"""Tests for Merzoug Assembly.

The method's promise is narrow and testable: annealing each oligo pair and
ligating them in order must reproduce the designed construct exactly, with
no PCR and no restriction digestion at internal junctions. Everything below
tests that promise rather than the implementation's shape.
"""
import random

import pytest

from gsynth_engine.constants import left_remainders, right_remainders
from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.sequence import (
    SequenceError,
    is_palindrome,
    reverse_complement,
)

# A realistic peptide-coding insert, long enough to need several fragments.
LONG_INSERT = (
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGGGCGGC"
    "GGCCCGGGCGCGGGCAGCCTGCAGCCGCTGGCGCTGGAAGGCAGCCTGCAGAAACGCGGCATCGTGGAA"
    "CAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACTAA"
)


def random_insert(length: int, seed: int) -> str:
    rng = random.Random(seed)
    # Keep it in frame and stop-free so it reads like a real coding insert.
    codons = [
        c for c in
        ("GCG", "TGC", "GAT", "GAA", "TTT", "GGC", "CAT", "ATT", "AAA", "CTG",
         "ATG", "AAC", "CCG", "CAG", "CGC", "AGC", "ACC", "GTG", "TGG", "TAT")
    ]
    return "".join(rng.choice(codons) for _ in range(length // 3)) + "TAA"


class TestAssemblyReproducesTheConstruct:
    """The one property that matters: ligation rebuilds the design."""

    def test_forward_oligos_rebuild_the_top_strand(self):
        plan = design_merzoug_assembly(LONG_INSERT)
        rebuilt = "".join(f.forward for f in plan.fragments)
        assert rebuilt == plan.construct_forward

    def test_reverse_oligos_rebuild_the_bottom_strand(self):
        """Each reverse oligo flipped back into top-strand sense, in fragment
        order, must rebuild the construct's bottom strand."""
        plan = design_merzoug_assembly(LONG_INSERT)
        rebuilt = "".join(reverse_complement(f.reverse) for f in plan.fragments)
        assert rebuilt == reverse_complement(plan.construct_reverse)

    def test_verify_reports_no_problems(self):
        plan = design_merzoug_assembly(LONG_INSERT)
        assert plan.verify() == []

    def test_construct_equals_the_ssd_design(self):
        """Fragmenting must not change the construct — only cut it up."""
        from gsynth_engine.ssd import design_small_sequence
        plan = design_merzoug_assembly(LONG_INSERT)
        ssd = design_small_sequence(LONG_INSERT)
        assert plan.construct_forward == ssd.forward
        assert plan.construct_reverse == ssd.reverse

    @pytest.mark.parametrize("length", [220, 400, 650, 900, 1500])
    @pytest.mark.parametrize("overhang", [4, 6, 8])
    def test_holds_across_lengths_and_overhangs(self, length, overhang):
        insert = random_insert(length, seed=length + overhang)
        plan = design_merzoug_assembly(
            insert, overhang_length=overhang, target_oligo_length=90,
        )
        assert plan.verify() == []
        assert "".join(f.forward for f in plan.fragments) == plan.construct_forward


class TestJunctions:
    """Overhang quality decides whether the assembly is ordered or scrambled."""

    def test_junction_overhangs_have_the_requested_length(self):
        for overhang in (4, 5, 6, 7, 8):
            plan = design_merzoug_assembly(LONG_INSERT, overhang_length=overhang)
            for junction in plan.junction_overhangs:
                assert len(junction) == overhang

    def test_junction_overhangs_are_unique(self):
        """Two identical overhangs let fragments ligate in the wrong order."""
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        junctions = plan.junction_overhangs
        assert len(junctions) == len(set(junctions))

    def test_junction_overhangs_are_not_palindromic(self):
        """A palindromic overhang anneals to itself — fragments self-ligate."""
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        for junction in plan.junction_overhangs:
            assert not is_palindrome(junction), junction

    def test_no_junction_reuses_a_terminal_enzyme_overhang(self):
        """Otherwise a fragment could ligate straight into the cut vector."""
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        terminal = {plan.ssd.left_overhang, plan.ssd.right_overhang}
        for junction in plan.junction_overhangs:
            assert junction not in terminal
            assert reverse_complement(junction) not in terminal

    def test_adjacent_fragments_share_a_complementary_junction(self):
        """Fragment i's right overhang must be fragment i+1's left overhang."""
        plan = design_merzoug_assembly(LONG_INSERT)
        for left, right in zip(plan.fragments, plan.fragments[1:]):
            assert left.right_overhang == right.left_overhang
            assert len(left.right_overhang) == plan.overhang_length

    def test_overhangs_are_not_all_at_or_all_gc(self):
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        for junction in plan.junction_overhangs:
            gc = sum(1 for base in junction if base in "GC")
            assert 0 < gc < len(junction), junction


class TestTerminalEnds:
    """The two ends must match the vector, not each other."""

    def test_first_fragment_carries_the_left_enzyme_overhang(self):
        plan = design_merzoug_assembly(LONG_INSERT, enzyme_pair="NdeI / XhoI")
        first = plan.fragments[0]
        assert first.is_first
        assert first.left_overhang == "TA"
        assert first.forward.startswith(left_remainders("NdeI")[0])   # TATG

    def test_last_fragment_carries_the_right_enzyme_overhang(self):
        plan = design_merzoug_assembly(LONG_INSERT, enzyme_pair="NdeI / XhoI")
        last = plan.fragments[-1]
        assert last.is_last
        assert last.right_overhang == "TCGA"
        assert last.reverse.startswith(right_remainders("XhoI")[1])   # TCGAG

    def test_works_with_another_enzyme_pair(self):
        plan = design_merzoug_assembly(
            LONG_INSERT, enzyme_pair="BamHI / EcoRI", cleavage_site="TEV",
        )
        assert plan.verify() == []
        assert plan.fragments[0].left_overhang == "GATC"     # BamHI
        assert plan.fragments[-1].right_overhang == "AATT"   # EcoRI


class TestOligoSizing:
    def test_oligos_land_near_the_requested_length(self):
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=90)
        for fragment in plan.fragments:
            assert 40 <= len(fragment.forward) <= 160, len(fragment.forward)

    def test_smaller_target_gives_more_fragments(self):
        few = design_merzoug_assembly(LONG_INSERT, target_oligo_length=150)
        many = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        assert many.fragment_count > few.fragment_count

    def test_short_insert_needs_no_fragmentation(self):
        plan = design_merzoug_assembly("GGCATCGTGGAACAGTGCTGCACCAGCTAA")
        assert plan.fragment_count == 1
        assert plan.fragments[0].is_first and plan.fragments[0].is_last
        assert plan.verify() == []

    def test_oligo_count_is_two_per_fragment(self):
        plan = design_merzoug_assembly(LONG_INSERT)
        assert plan.oligo_count == 2 * plan.fragment_count


class TestInputErrors:
    """Errors must say what to do, not just what went wrong."""

    def test_overhang_outside_the_method_is_refused(self):
        with pytest.raises(SequenceError, match="between 4 and 8"):
            design_merzoug_assembly(LONG_INSERT, overhang_length=2)
        with pytest.raises(SequenceError, match="between 4 and 8"):
            design_merzoug_assembly(LONG_INSERT, overhang_length=12)

    def test_target_shorter_than_the_overhangs_is_refused(self):
        with pytest.raises(SequenceError, match="too short"):
            design_merzoug_assembly(LONG_INSERT, target_oligo_length=10, overhang_length=6)

    def test_invalid_bases_are_reported(self):
        with pytest.raises(SequenceError, match="not A, C, G or T"):
            design_merzoug_assembly("ATGXYZATG")

    def test_empty_input_is_reported(self):
        with pytest.raises(SequenceError, match="empty"):
            design_merzoug_assembly("   ")
