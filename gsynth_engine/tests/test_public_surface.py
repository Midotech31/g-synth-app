"""The values the interface reads, and the sentences the user is shown.

Two things that look like plumbing and are not.

The properties on `SSDResult` and `OligoPair` go straight into the API
payload — `forward_tm` is the temperature printed beside an oligo on the
order sheet, `right_overhang_strand` decides which way a duplex is drawn.
Nothing asserted any of them; they were exercised only through whatever a
higher-level test happened to look at.

And `SequenceError` messages are shown verbatim in the interface. A message
that names the problem without naming the fix costs the user a round trip.
"""
import pytest

from gsynth_engine.constants import HIS_TAG, LEFT_LINKER, RIGHT_LINKER
from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.thermo import ANNEALING, melting_temperature

INSERT = ("GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGC"
          "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCTAA")


class TestTheCassetteCombinations:
    """Four combinations of tag and linkers; three had no test.

    Each is a checkbox pair in the interface, and each changes the protein.
    """

    def test_tag_and_linkers_together(self):
        r = design_small_sequence(INSERT, include_his_tag=True, include_linkers=True)
        assert LEFT_LINKER + HIS_TAG + RIGHT_LINKER in r.forward
        assert "6×His tag" in [s.name for s in r.segments]

    def test_tag_without_linkers(self):
        """A tag butted straight onto the protein — no flexible spacer. The
        user asked for exactly that; it must not silently add one."""
        r = design_small_sequence(INSERT, include_his_tag=True, include_linkers=False)
        assert HIS_TAG in r.forward
        assert LEFT_LINKER not in r.forward
        assert RIGHT_LINKER not in r.forward

    def test_linkers_without_a_tag(self):
        """A spacer with nothing to space. Unusual, but it is what the two
        boxes allow, and it must not quietly insert a tag."""
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
        """The reverse strand is built by a parallel set of branches, so a
        combination can be right on top and wrong underneath."""
        r = design_small_sequence(INSERT, include_his_tag=tag, include_linkers=linkers)
        plan = design_merzoug_assembly(INSERT, include_his_tag=tag,
                                       include_linkers=linkers)
        assert plan.verify() == []
        assert r.forward == plan.construct_forward


class TestTheNumbersTheOrderSheetPrints:
    """These are read off a screen and typed into a supplier's form."""

    def test_ssd_lengths_and_gc_match_the_strands_themselves(self):
        r = design_small_sequence(INSERT)
        assert r.forward_length == len(r.forward)
        assert r.reverse_length == len(r.reverse)
        assert r.forward_gc == round(gc_content(r.forward), 1)
        assert r.reverse_gc == round(gc_content(r.reverse), 1)

    def test_ssd_tm_is_the_annealing_reaction_not_a_generic_dilution(self):
        """The two differ by about 7 °C, and the protocol anneals at the
        first. A Tm quoted at the wrong concentration is a failed anneal."""
        r = design_small_sequence(INSERT)
        assert r.forward_tm == round(
            melting_temperature(r.forward, conditions=ANNEALING), 1)
        assert r.reverse_tm == round(
            melting_temperature(r.reverse, conditions=ANNEALING), 1)

    def test_fragment_numbers_match_their_own_oligos(self):
        plan = design_merzoug_assembly(INSERT, target_oligo_length=60)
        for f in plan.fragments:
            assert f.forward_length == len(f.forward)
            assert f.reverse_length == len(f.reverse)
            assert f.forward_tm == round(
                melting_temperature(f.forward, conditions=ANNEALING), 1)
            assert f.duplex_gc == round(gc_content(f.forward), 1)


class TestWhichStrandCarriesEachOverhang:
    """The drawing is laid out from these. Getting one wrong puts an overhang
    on the strand that cannot present it, and the picture still looks fine."""

    def test_a_five_prime_pair_puts_them_on_opposite_strands(self):
        plan = design_merzoug_assembly(INSERT, enzyme_pair="NdeI / XhoI")
        assert plan.fragments[0].left_overhang_strand == "top"
        assert plan.fragments[-1].right_overhang_strand == "bottom"

    def test_a_three_prime_pair_is_the_other_way_round(self):
        """KpnI and SacI leave 3' overhangs — the same protrusion sits on the
        other strand, which is the distinction the whole model turns on."""
        plan = design_merzoug_assembly(INSERT, enzyme_pair="KpnI / SacI")
        assert plan.fragments[0].left_overhang_strand == "bottom"
        assert plan.fragments[-1].right_overhang_strand == "top"

    def test_a_blunt_end_says_blunt(self):
        plan = design_merzoug_assembly(INSERT, enzyme_pair="EcoRV / SmaI")
        assert plan.fragments[0].left_overhang_strand == "blunt"
        assert plan.fragments[-1].right_overhang_strand == "blunt"

    def test_internal_junctions_are_always_five_prime_on_top(self):
        """However the outer ends are cut, the stagger between fragments is
        made by the design and is always the same way round."""
        plan = design_merzoug_assembly(INSERT, enzyme_pair="KpnI / SacI",
                                       target_oligo_length=60)
        for fragment in plan.fragments[1:]:
            assert fragment.left_overhang_strand == "top"


class TestErrorsTellTheUserWhatToDo:
    """`SequenceError` text reaches the interface unchanged."""

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
