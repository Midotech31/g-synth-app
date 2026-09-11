import pytest

from gsynth_engine.constants import HIS_TAG, LEFT_LINKER, RIGHT_LINKER
from gsynth_engine.esd import design_extended_sequence
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.thermo import ANNEALING, melting_temperature

INSERT = ("GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
          "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCTAA")


class TestTheCassetteCombinations:


    def test_tag_and_linkers_together(self):
        r = design_small_sequence(INSERT, include_his_tag=True, include_linkers=True)
        assert LEFT_LINKER + HIS_TAG + RIGHT_LINKER in r.forward
        assert "6×His tag" in [s.name for s in r.segments]

    def test_tag_without_linkers(self):

        r = design_small_sequence(INSERT, include_his_tag=True, include_linkers=False)
        assert HIS_TAG in r.forward
        assert LEFT_LINKER not in r.forward
        assert RIGHT_LINKER not in r.forward

    def test_linkers_without_a_tag(self):

        r = design_small_sequence(INSERT, include_his_tag=False, include_linkers=True)
        assert LEFT_LINKER + RIGHT_LINKER in r.forward
        assert HIS_TAG not in r.forward

    def test_neither(self):
        r = design_small_sequence(INSERT, include_his_tag=False, include_linkers=False)
        assert HIS_TAG not in r.forward
        assert LEFT_LINKER not in r.forward

    @pytest.mark.parametrize("tag,linkers", [(True, True), (True, False),
                                             (False, True), (False, False)])
    def test_both_strands_stay_complementary_in_every_combination(self, tag, linkers):

        r = design_small_sequence(INSERT, include_his_tag=tag, include_linkers=linkers)
        plan = design_extended_sequence(INSERT, include_his_tag=tag,
                                       include_linkers=linkers)
        assert plan.verify() == []
        assert r.forward == plan.construct_forward


class TestTheNumbersTheOrderSheetPrints:


    def test_ssd_lengths_and_gc_match_the_strands_themselves(self):
        r = design_small_sequence(INSERT)
        assert r.forward_length == len(r.forward)
        assert r.reverse_length == len(r.reverse)
        assert r.forward_gc == round(gc_content(r.forward), 1)
        assert r.reverse_gc == round(gc_content(r.reverse), 1)

    def test_ssd_tm_is_the_annealing_reaction_not_a_generic_dilution(self):

        r = design_small_sequence(INSERT)
        assert r.forward_tm == round(
            melting_temperature(r.forward, conditions=ANNEALING), 1)
        assert r.reverse_tm == round(
            melting_temperature(r.reverse, conditions=ANNEALING), 1)

    def test_fragment_numbers_match_their_own_oligos(self):
        plan = design_extended_sequence(INSERT, target_oligo_length=60)
        for f in plan.fragments:
            assert f.forward_length == len(f.forward)
            assert f.reverse_length == len(f.reverse)
            assert f.forward_tm == round(
                melting_temperature(f.forward, conditions=ANNEALING), 1)
            assert f.duplex_gc == round(gc_content(f.forward), 1)


class TestWhichStrandCarriesEachOverhang:


    def test_a_five_prime_pair_puts_them_on_opposite_strands(self):
        plan = design_extended_sequence(INSERT, enzyme_pair="NdeI / XhoI")
        assert plan.fragments[0].left_overhang_strand == "top"
        assert plan.fragments[-1].right_overhang_strand == "bottom"

    def test_a_three_prime_pair_is_the_other_way_round(self):

        plan = design_extended_sequence(INSERT, enzyme_pair="KpnI / SacI")
        assert plan.fragments[0].left_overhang_strand == "bottom"
        assert plan.fragments[-1].right_overhang_strand == "top"

    def test_a_blunt_end_says_blunt(self):
        plan = design_extended_sequence(INSERT, enzyme_pair="EcoRV / SmaI")
        assert plan.fragments[0].left_overhang_strand == "blunt"
        assert plan.fragments[-1].right_overhang_strand == "blunt"

    def test_internal_junctions_are_always_five_prime_on_top(self):

        plan = design_extended_sequence(INSERT, enzyme_pair="KpnI / SacI",
                                       target_oligo_length=60)
        for fragment in plan.fragments[1:]:
            assert fragment.left_overhang_strand == "top"


class TestErrorsTellTheUserWhatToDo:


    def test_an_unknown_enzyme_lists_the_ones_that_are_known(self):
        with pytest.raises(SequenceError) as caught:
            design_small_sequence(INSERT, enzyme_pair="BsaI / XhoI")
        message = str(caught.value)
        assert "BsaI" in message
        assert "NdeI" in message, "the message must name enzymes that do work"
        assert "more" in message, "and say how many others exist"

    def test_an_unknown_protease_site_lists_the_alternatives(self):
        with pytest.raises(SequenceError) as caught:
            design_small_sequence(INSERT, cleavage_site="Papain")
        message = str(caught.value)
        assert "Thrombin" in message and "TEV" in message

    def test_a_malformed_enzyme_pair_says_what_the_format_is(self):
        with pytest.raises(SequenceError) as caught:
            design_small_sequence(INSERT, enzyme_pair="NdeI")
        assert "/" in str(caught.value)
