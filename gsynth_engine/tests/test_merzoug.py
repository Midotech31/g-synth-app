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

    def test_no_two_junctions_differ_by_a_single_base(self):
        """T4 ligase cross-ligates near-identical overhangs.

        NEB's ligase-fidelity data show that overhangs differing at one
        position join at a measurable rate. Two such junctions in the same
        reaction produce a misassembled construct that looks fine on a gel
        and only shows up at sequencing, so they are kept apart by design.
        """
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        junctions = plan.junction_overhangs
        assert len(junctions) >= 2, "need several junctions for this to mean anything"

        for i, first in enumerate(junctions):
            for second in junctions[i + 1 :]:
                for candidate in (second, reverse_complement(second)):
                    mismatches = sum(1 for a, b in zip(first, candidate) if a != b)
                    assert mismatches >= 2, (
                        f"{first} and {second} differ by {mismatches} base(s)"
                    )

    def test_fidelity_rule_also_covers_the_terminal_overhangs(self):
        """A junction one base from the vector's sticky end can ligate into it."""
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=60)
        terminal = [plan.ssd.left_overhang, plan.ssd.right_overhang]

        for junction in plan.junction_overhangs:
            for end in terminal:
                if len(end) != len(junction):
                    continue
                for candidate in (end, reverse_complement(end)):
                    mismatches = sum(1 for a, b in zip(junction, candidate) if a != b)
                    assert mismatches >= 2, f"{junction} too close to terminal {end}"


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


class TestLongConstructs:
    """The method exists to build genes, not peptides.

    Every test above uses an insert of a few hundred bases, and for a long
    time that was all the design was ever run on. A 2.4 kb gene — an ordinary
    target — cut into 90 nt oligos needs 26 junctions, and there are only 22
    sets of 4 nt overhangs that ligate in one order no matter how they are
    chosen. Designing it used to fail outright at junction 19.
    """

    def test_the_supply_table_matches_the_rules_it_describes(self):
        """OVERHANG_SUPPLY is measured, and the measurement can go stale.

        The table is what decides when to widen. Tighten the rules in
        `_OverhangPool` without re-deriving it and the design will keep
        promising a supply of overhangs that no longer exists, then fail at
        the far end of a long gene — exactly the failure it is there to
        prevent.
        """
        import itertools

        from gsynth_engine.merzoug import OVERHANG_SUPPLY, _OverhangPool

        for length, promised in OVERHANG_SUPPLY.items():
            pool = _OverhangPool(set())
            for word in itertools.product("ACGT", repeat=length):
                candidate = "".join(word)
                if pool.problem(candidate) is None:
                    pool.take(candidate)
            assert len(pool.taken) == promised, (
                f"{length} nt overhangs now supply {len(pool.taken)} junctions, "
                f"not the {promised} the table promises."
            )

    @pytest.mark.parametrize("length", [2400, 5000])
    def test_a_whole_gene_designs_and_re_ligates(self, length):
        """The load-bearing property must hold at gene scale, not just peptide
        scale: the fragments still rebuild the construct on both strands."""
        plan = design_merzoug_assembly(random_insert(length, seed=7),
                                       target_oligo_length=90)
        assert plan.fragment_count > 25
        assert plan.verify() == []
        assert "".join(f.forward for f in plan.fragments) == plan.construct_forward
        assert "".join(
            reverse_complement(f.reverse) for f in plan.fragments
        ) == reverse_complement(plan.construct_reverse)

    def test_overhangs_widen_when_4_nt_cannot_supply_the_junctions(self):
        plan = design_merzoug_assembly(random_insert(2400, seed=7),
                                       target_oligo_length=90)
        assert plan.overhang_length > 4
        assert all(
            len(o) == plan.overhang_length for o in plan.junction_overhangs
        )

    def test_the_widening_is_reported_rather_than_silent(self):
        """The user asked for 4 nt overhangs and will get something else.

        Ordering oligos built on an assumption the design quietly abandoned is
        how a bench protocol stops matching the tubes on the bench.
        """
        plan = design_merzoug_assembly(random_insert(2400, seed=7),
                                       target_oligo_length=90)
        assert any("widened to" in w for w in plan.warnings)

    def test_widening_does_not_lengthen_the_oligos(self):
        """Fragments are cut from a fixed construct, so a wider stagger moves
        where the cuts fall, not how much has to be synthesised. If this ever
        stops being true, widening starts costing money per junction."""
        plan = design_merzoug_assembly(random_insert(2400, seed=7),
                                       target_oligo_length=90)
        assert plan.longest_oligo <= 130

    def test_a_short_insert_still_gets_the_4_nt_overhangs_asked_for(self):
        """Widening must be a last resort. A design that fits in 4 nt has to
        keep them, or every existing protocol in the lab changes at once."""
        plan = design_merzoug_assembly(LONG_INSERT, target_oligo_length=90)
        assert plan.overhang_length == 4
        assert not any("widened" in w for w in plan.warnings)

    def test_junctions_stay_mutually_distinct_at_gene_scale(self):
        """Uniqueness is what fixes the assembly order. Widening enlarges the
        supply; it must not relax the rule the supply exists to satisfy."""
        plan = design_merzoug_assembly(random_insert(5000, seed=11),
                                       target_oligo_length=90)
        overhangs = plan.junction_overhangs
        assert len(set(overhangs)) == len(overhangs)
        for overhang in overhangs:
            assert not is_palindrome(overhang)
            assert reverse_complement(overhang) not in set(overhangs) - {overhang}

    def test_a_tandem_repeat_is_refused_by_name(self):
        """A repeat cannot be assembled by overhang-directed ligation at all —
        the whole gene contains only as many distinct words as its unit is
        long. Blaming the search sends the user to change settings that cannot
        help; the message has to name the sequence.
        """
        repeat = "ATG" + "GGTCCGGCTGGTCCGGCT" * 50 + "TAA"
        with pytest.raises(SequenceError, match="not contain enough distinct"):
            design_merzoug_assembly(repeat, target_oligo_length=90)

    def test_a_long_ordinary_gene_is_not_called_repetitive(self):
        """The measure of "too few distinct words" is what the junctions
        consume, not the gene's length.

        No sequence of any kind has more than 4^8 distinct 8 nt words, so a
        rule phrased against length calls every gene past 131 kb repetitive —
        including one drawn base by base at random.
        """
        import random as _random

        rng = _random.Random(17)
        ordinary = "".join(rng.choice("ACGT") for _ in range(150_000))
        plan = design_merzoug_assembly(ordinary, target_oligo_length=90)
        assert plan.verify() == []

    def test_placement_does_not_grow_quadratically(self):
        """Junction placement was O(junctions²) and the endpoint accepts 200 kb.

        Every candidate overhang was compared against every one already
        placed, which is invisible on a peptide and takes about half a minute
        on the largest input the API will accept — a hang any signed-in user
        could trigger by accident. The ceiling here is thirty times the
        measured cost, so it fails for a return to quadratic and not for a
        slow machine.
        """
        import random as _random
        import time

        # Maximal word diversity, so this measures placement rather than the
        # sequence running out of distinct overhangs.
        rng = _random.Random(3)
        longest = "".join(rng.choice("ACGT") for _ in range(200_000))
        started = time.perf_counter()
        plan = design_merzoug_assembly(longest, target_oligo_length=90)
        elapsed = time.perf_counter() - started

        assert plan.verify() == []
        assert elapsed < 5.0, f"placement took {elapsed:.1f}s"
