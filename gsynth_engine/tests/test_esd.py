import random

import pytest

from gsynth_engine.constants import (
    RESTRICTION_ENZYMES,
    left_remainders,
    right_remainders,
)
from gsynth_engine.constants import overhang as enzyme_overhang
from gsynth_engine.esd import design_extended_sequence
from gsynth_engine.sequence import (
    SequenceError,
    is_palindrome,
    reverse_complement,
)

LONG_INSERT = (
    "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
    "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
    "TTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGGGCGGC"
    "GGCCCGGGCGCGGGCAGCCTGCAGCCGCTGGCGCTGGAAGGCAGCCTGCAGAAACGCGGCATCGTGGAA"
    "CAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACTAA"
)


def random_insert(length: int, seed: int) -> str:
    rng = random.Random(seed)

    codons = [
        c for c in
        ("GCG", "TGC", "GAT", "GAA", "TTT", "GGC", "CAT", "ATT", "AAA", "CTG",
         "ATG", "AAC", "CCG", "CAG", "CGC", "AGC", "ACC", "GTG", "TGG", "TAT")
    ]
    return "".join(rng.choice(codons) for _ in range(length // 3)) + "TAA"


class TestAssemblyReproducesTheConstruct:


    def test_forward_oligos_rebuild_the_top_strand(self):
        plan = design_extended_sequence(LONG_INSERT)
        rebuilt = "".join(f.forward for f in plan.fragments)
        assert rebuilt == plan.construct_forward

    def test_reverse_oligos_rebuild_the_bottom_strand(self):

        plan = design_extended_sequence(LONG_INSERT)
        rebuilt = "".join(reverse_complement(f.reverse) for f in plan.fragments)
        assert rebuilt == reverse_complement(plan.construct_reverse)

    def test_verify_reports_no_problems(self):
        plan = design_extended_sequence(LONG_INSERT)
        assert plan.verify() == []

    def test_construct_equals_the_ssd_design(self):

        from gsynth_engine.ssd import design_small_sequence
        plan = design_extended_sequence(LONG_INSERT)
        ssd = design_small_sequence(LONG_INSERT)
        assert plan.construct_forward == ssd.forward
        assert plan.construct_reverse == ssd.reverse

    @pytest.mark.parametrize("length", [220, 400, 650, 900, 1500])
    @pytest.mark.parametrize("overhang", [4, 6, 8])
    def test_holds_across_lengths_and_overhangs(self, length, overhang):
        insert = random_insert(length, seed=length + overhang)
        plan = design_extended_sequence(
            insert, overhang_length=overhang, target_oligo_length=90,
        )
        assert plan.verify() == []
        assert "".join(f.forward for f in plan.fragments) == plan.construct_forward


class TestJunctions:


    def test_junction_overhangs_have_the_requested_length(self):
        for overhang in (4, 5, 6, 7, 8):
            plan = design_extended_sequence(LONG_INSERT, overhang_length=overhang)
            for junction in plan.junction_overhangs:
                assert len(junction) == overhang

    def test_junction_overhangs_are_unique(self):

        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        junctions = plan.junction_overhangs
        assert len(junctions) == len(set(junctions))

    def test_junction_overhangs_are_not_palindromic(self):

        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        for junction in plan.junction_overhangs:
            assert not is_palindrome(junction), junction

    def test_no_junction_reuses_a_terminal_enzyme_overhang(self):

        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        terminal = {plan.ssd.left_overhang, plan.ssd.right_overhang}
        for junction in plan.junction_overhangs:
            assert junction not in terminal
            assert reverse_complement(junction) not in terminal

    def test_adjacent_fragments_share_a_complementary_junction(self):

        plan = design_extended_sequence(LONG_INSERT)
        for left, right in zip(plan.fragments, plan.fragments[1:], strict=False):
            assert left.right_overhang == right.left_overhang
            assert len(left.right_overhang) == plan.overhang_length

    def test_overhangs_are_not_all_at_or_all_gc(self):
        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        for junction in plan.junction_overhangs:
            gc = sum(1 for base in junction if base in "GC")
            assert 0 < gc < len(junction), junction

    def test_no_two_junctions_differ_by_a_single_base(self):

        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        junctions = plan.junction_overhangs
        assert len(junctions) >= 2, "need several junctions for this to mean anything"

        for i, first in enumerate(junctions):
            for second in junctions[i + 1 :]:
                for candidate in (second, reverse_complement(second)):
                    mismatches = sum(1 for a, b in zip(first, candidate, strict=False) if a != b)
                    assert mismatches >= 2, (
                        f"{first} and {second} differ by {mismatches} base(s)"
                    )

    def test_fidelity_rule_also_covers_the_terminal_overhangs(self):

        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        terminal = [plan.ssd.left_overhang, plan.ssd.right_overhang]

        for junction in plan.junction_overhangs:
            for end in terminal:
                if len(end) != len(junction):
                    continue
                for candidate in (end, reverse_complement(end)):
                    mismatches = sum(1 for a, b in zip(junction, candidate, strict=False) if a != b)
                    assert mismatches >= 2, f"{junction} too close to terminal {end}"


class TestTerminalEnds:


    def test_first_fragment_carries_the_left_enzyme_overhang(self):
        plan = design_extended_sequence(LONG_INSERT, enzyme_pair="NdeI / XhoI")
        first = plan.fragments[0]
        assert first.is_first
        assert first.left_overhang == "TA"
        assert first.forward.startswith(left_remainders("NdeI")[0])

    def test_last_fragment_carries_the_right_enzyme_overhang(self):
        plan = design_extended_sequence(LONG_INSERT, enzyme_pair="NdeI / XhoI")
        last = plan.fragments[-1]
        assert last.is_last
        assert last.right_overhang == "TCGA"
        assert last.reverse.startswith(right_remainders("XhoI")[1])

    def test_works_with_another_enzyme_pair(self):
        plan = design_extended_sequence(
            LONG_INSERT, enzyme_pair="BamHI / EcoRI", cleavage_site="TEV",
        )
        assert plan.verify() == []
        assert plan.fragments[0].left_overhang == "GATC"
        assert plan.fragments[-1].right_overhang == "AATT"


class TestOligoSizing:
    def test_oligos_land_near_the_requested_length(self):
        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=90)
        for fragment in plan.fragments:
            assert 40 <= len(fragment.forward) <= 160, len(fragment.forward)

    def test_smaller_target_gives_more_fragments(self):
        few = design_extended_sequence(LONG_INSERT, target_oligo_length=150)
        many = design_extended_sequence(LONG_INSERT, target_oligo_length=60)
        assert many.fragment_count > few.fragment_count

    def test_short_insert_needs_no_fragmentation(self):
        plan = design_extended_sequence("GGCATCGTGGAACAGTGCTGCACCAGCTAA")
        assert plan.fragment_count == 1
        assert plan.fragments[0].is_first and plan.fragments[0].is_last
        assert plan.verify() == []

    def test_oligo_count_is_two_per_fragment(self):
        plan = design_extended_sequence(LONG_INSERT)
        assert plan.oligo_count == 2 * plan.fragment_count


class TestInputErrors:


    def test_overhang_outside_the_method_is_refused(self):
        with pytest.raises(SequenceError, match="between 4 and 8"):
            design_extended_sequence(LONG_INSERT, overhang_length=2)
        with pytest.raises(SequenceError, match="between 4 and 8"):
            design_extended_sequence(LONG_INSERT, overhang_length=12)

    def test_target_shorter_than_the_overhangs_is_refused(self):
        with pytest.raises(SequenceError, match="too short"):
            design_extended_sequence(LONG_INSERT, target_oligo_length=10, overhang_length=6)

    def test_invalid_bases_are_reported(self):
        with pytest.raises(SequenceError, match="not A, C, G or T"):
            design_extended_sequence("ATGXYZATG")

    def test_empty_input_is_reported(self):
        with pytest.raises(SequenceError, match="empty"):
            design_extended_sequence("   ")


class TestLongConstructs:


    def test_the_supply_table_matches_the_rules_it_describes(self):

        import itertools

        from gsynth_engine.esd import OVERHANG_SUPPLY, _OverhangPool

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

        plan = design_extended_sequence(random_insert(length, seed=7),
                                       target_oligo_length=90)
        assert plan.fragment_count > 25
        assert plan.verify() == []
        assert "".join(f.forward for f in plan.fragments) == plan.construct_forward
        assert "".join(
            reverse_complement(f.reverse) for f in plan.fragments
        ) == reverse_complement(plan.construct_reverse)

    def test_overhangs_widen_when_4_nt_cannot_supply_the_junctions(self):
        plan = design_extended_sequence(random_insert(2400, seed=7),
                                       target_oligo_length=90)
        assert plan.overhang_length > 4
        assert all(
            len(o) == plan.overhang_length for o in plan.junction_overhangs
        )

    def test_the_widening_is_reported_rather_than_silent(self):

        plan = design_extended_sequence(random_insert(2400, seed=7),
                                       target_oligo_length=90)
        assert any("widened to" in w for w in plan.warnings)

    def test_widening_does_not_lengthen_the_oligos(self):

        plan = design_extended_sequence(random_insert(2400, seed=7),
                                       target_oligo_length=90)
        assert plan.longest_oligo <= 130

    def test_a_short_insert_still_gets_the_4_nt_overhangs_asked_for(self):

        plan = design_extended_sequence(LONG_INSERT, target_oligo_length=90)
        assert plan.overhang_length == 4
        assert not any("widened" in w for w in plan.warnings)

    def test_junctions_stay_mutually_distinct_at_gene_scale(self):

        plan = design_extended_sequence(random_insert(5000, seed=11),
                                       target_oligo_length=90)
        overhangs = plan.junction_overhangs
        assert len(set(overhangs)) == len(overhangs)
        for overhang in overhangs:
            assert not is_palindrome(overhang)
            assert reverse_complement(overhang) not in set(overhangs) - {overhang}

    def test_a_tandem_repeat_is_refused_by_name(self):

        repeat = "ATG" + "GGTCCGGCTGGTCCGGCT" * 50 + "TAA"
        with pytest.raises(SequenceError, match="not contain enough distinct"):
            design_extended_sequence(repeat, target_oligo_length=90)

    def test_a_long_ordinary_gene_is_not_called_repetitive(self):

        import random as _random

        rng = _random.Random(17)
        ordinary = "".join(rng.choice("ACGT") for _ in range(150_000))
        plan = design_extended_sequence(ordinary, target_oligo_length=90)
        assert plan.verify() == []

    def test_placement_does_not_grow_quadratically(self):

        import random as _random
        import time


        rng = _random.Random(3)
        longest = "".join(rng.choice("ACGT") for _ in range(200_000))
        started = time.perf_counter()
        plan = design_extended_sequence(longest, target_oligo_length=90)
        elapsed = time.perf_counter() - started

        assert plan.verify() == []
        assert elapsed < 5.0, f"placement took {elapsed:.1f}s"


class TestTheAssembledEndsMatchTheEnzymes:


    def _observed(self, plan):

        return plan.terminal_ends

    @pytest.mark.parametrize("enzyme", sorted(RESTRICTION_ENZYMES))
    def test_left_end_is_what_the_left_enzyme_leaves(self, enzyme):
        plan = design_extended_sequence(
            LONG_INSERT, enzyme_pair=f"{enzyme} / XhoI", target_oligo_length=90,
        )
        left, _ = self._observed(plan)
        assert left == enzyme_overhang(enzyme)
        assert plan.ssd.left_overhang == enzyme_overhang(enzyme)[0]

    @pytest.mark.parametrize("enzyme", sorted(RESTRICTION_ENZYMES))
    def test_right_end_is_what_the_right_enzyme_leaves(self, enzyme):
        plan = design_extended_sequence(
            LONG_INSERT, enzyme_pair=f"NdeI / {enzyme}", target_oligo_length=90,
        )
        _, right = self._observed(plan)
        assert right == enzyme_overhang(enzyme)
        assert plan.ssd.right_overhang == enzyme_overhang(enzyme)[0]

    def test_a_three_prime_cutter_is_reported_as_one(self):

        plan = design_extended_sequence(
            LONG_INSERT, enzyme_pair="KpnI / SacI", target_oligo_length=90,
        )
        left, right = self._observed(plan)
        assert left == ("GTAC", "3'")
        assert right == ("AGCT", "3'")

    def test_ends_survive_a_gene_sized_assembly(self):

        plan = design_extended_sequence(
            random_insert(2_400, seed=7), target_oligo_length=90,
        )
        assert plan.overhang_length > 4, "this case should have widened"
        assert plan.terminal_ends == (enzyme_overhang("NdeI"), enzyme_overhang("XhoI"))
        assert plan.verify() == []

    def test_verify_catches_an_end_that_does_not_match(self):

        from dataclasses import replace

        plan = design_extended_sequence(LONG_INSERT, enzyme_pair="NdeI / XhoI")
        assert plan.verify() == []

        first = plan.fragments[0]
        plan.fragments[0] = replace(first, bottom_offset=first.bottom_offset + 1)

        problems = plan.verify()
        assert any("left-hand end" in p for p in problems), problems
        assert any("NdeI" in p for p in problems), problems
