from __future__ import annotations

import pytest

from gsynth_engine.constants import RESTRICTION_ENZYMES
from gsynth_engine.duplex import (
    GAP,
    construct_duplex,
    fragment_duplex,
    junction_view,
)
from gsynth_engine.esd import design_extended_sequence
from gsynth_engine.sequence import SequenceError, complement, reverse_complement

INSERT = "ATG" + "GCTAGCAAAGGTTTCCGTGAAGATCTGGCAAAATTCCTGCAGGCTAACGGT" * 3


@pytest.fixture(scope="module")
def plan():
    return design_extended_sequence(INSERT, target_oligo_length=70)


@pytest.fixture(scope="module")
def view(plan):
    return construct_duplex(plan)


class TestPairing:


    def test_no_mismatches_anywhere(self, view):
        assert view.mismatches() == []

    def test_every_paired_column_is_complementary(self, view):
        for top, bottom in zip(view.top, view.bottom, strict=False):
            if top != GAP and bottom != GAP:
                assert complement(top) == bottom

    def test_pairing_marks_match_the_strands(self, view):
        marks = view.paired()
        assert len(marks) == view.width
        for mark, top, bottom in zip(marks, view.top, view.bottom, strict=False):
            if top == GAP or bottom == GAP:
                assert mark == GAP
            else:
                assert mark == "|"

    def test_holds_for_every_enzyme_pair(self):

        names = sorted(RESTRICTION_ENZYMES)
        for left, right in zip(names, names[1:], strict=False):
            if left == right:
                continue
            design = design_extended_sequence(
                INSERT, enzyme_pair=f"{left} / {right}", target_oligo_length=70
            )
            assert construct_duplex(design).mismatches() == [], f"{left} / {right}"


class TestOverhangs:


    def test_left_overhang_is_top_strand_only(self, view, plan):

        width = len(plan.ssd.left_overhang)
        assert view.top[:width] == plan.ssd.left_overhang
        assert view.bottom[:width] == GAP * width
        assert view.bottom[width] != GAP

    def test_right_overhang_is_bottom_strand_only(self, view, plan):
        width = len(plan.ssd.right_overhang)
        assert view.top[-width:] == GAP * width
        assert view.top[-width - 1] != GAP

        assert view.bottom[-width:] == complement(plan.ssd.right_overhang)

    def test_a_three_prime_overhang_puts_the_stagger_on_the_other_strand(self):

        design = design_extended_sequence(
            INSERT, enzyme_pair="ApaI / XhoI", target_oligo_length=70
        )
        three_prime = construct_duplex(design)

        assert three_prime.mismatches() == []
        assert three_prime.bottom[0] != GAP
        assert three_prime.top[0] == GAP
        assert design.fragments[0].left_overhang_strand == "bottom"

    def test_frame_spans_both_strands_exactly(self, view, plan):
        bottom_sense = reverse_complement(plan.construct_reverse)
        assert view.top.strip(GAP) == plan.construct_forward
        assert view.bottom.strip(GAP) == complement(bottom_sense)

        assert not any(t == GAP and b == GAP for t, b in zip(view.top, view.bottom, strict=False))

    def test_blunt_ends_leave_no_stagger(self):

        design = design_extended_sequence(
            INSERT, enzyme_pair="EcoRV / EcoRV", target_oligo_length=70
        )
        blunt = construct_duplex(design)
        assert blunt.top[0] != GAP and blunt.bottom[0] != GAP
        assert blunt.top[-1] != GAP and blunt.bottom[-1] != GAP
        assert blunt.mismatches() == []


class TestFragmentSpans:


    def test_top_spans_tile_the_top_strand(self, view, plan):
        spans = view.top_fragments
        assert spans[0].start == 0
        assert spans[-1].end == len(plan.construct_forward)
        for left, right in zip(spans, spans[1:], strict=False):
            assert left.end == right.start

    def test_bottom_spans_tile_the_bottom_strand(self, view, plan):
        spans = view.bottom_fragments
        assert spans[0].start == len(plan.ssd.left_overhang)
        assert spans[-1].end == view.width
        for left, right in zip(spans, spans[1:], strict=False):
            assert left.end == right.start

    def test_the_two_strands_are_staggered_at_every_junction(self, view, plan):

        top_cuts = {span.end for span in view.top_fragments[:-1]}
        bottom_cuts = {span.end for span in view.bottom_fragments[:-1]}
        assert top_cuts.isdisjoint(bottom_cuts)

        for cut in sorted(top_cuts):
            offset = min(b - cut for b in bottom_cuts if b > cut)
            assert offset == plan.overhang_length

    def test_span_sequences_are_the_oligos_that_get_ordered(self, view, plan):
        for span, fragment in zip(view.top_fragments, plan.fragments, strict=False):
            assert view.top[span.start : span.end] == fragment.forward
        for span, fragment in zip(view.bottom_fragments, plan.fragments, strict=False):
            drawn = view.bottom[span.start : span.end]
            assert drawn == complement(reverse_complement(fragment.reverse))

    def test_junction_coordinates_are_the_top_strand_cuts(self, view, plan):
        assert view.junctions == [f.top_end for f in plan.fragments[:-1]]
        assert len(view.junctions) == plan.fragment_count - 1


class TestCassetteSegments:


    def test_segments_index_the_top_strand(self, view, plan):
        assert view.segments, "the cassette should be labelled"
        for span, segment in zip(view.segments, plan.ssd.segments, strict=False):
            assert span.name == segment.name
            assert view.top[span.start : span.end] == segment.sequence

    def test_segments_tile_the_top_strand_without_gaps(self, view, plan):
        assert view.segments[0].start == 0
        assert view.segments[-1].end == len(plan.construct_forward)
        for left, right in zip(view.segments, view.segments[1:], strict=False):
            assert left.end == right.start


class TestRows:


    def test_rows_reassemble_the_frame(self, view):
        rows = view.rows(60)
        assert "".join(row.top for row in rows) == view.top
        assert "".join(row.bottom for row in rows) == view.bottom

    def test_row_starts_advance_by_the_width(self, view):
        rows = view.rows(50)
        assert [row.start for row in rows] == list(range(0, view.width, 50))

    def test_numbering_counts_bases_not_columns(self, view):

        rows = view.rows(60)
        first = rows[0]
        assert first.top_start == 1
        assert first.bottom_start == 1

        second = rows[1]
        assert second.top_start == 61
        assert second.bottom_start == 61 - len(view.left_overhang)

    def test_row_end_numbers_are_the_last_base_in_the_row(self, view):
        for row in view.rows(60):
            if row.top_end is not None:
                assert row.top_end == row.top_start + sum(
                    1 for ch in row.top if ch != GAP
                ) - 1

    def test_a_row_with_no_bases_on_a_strand_reports_none(self, view):

        rows = view.rows(4)
        last = rows[-1]
        assert last.top_start is None
        assert last.bottom_start is not None


class TestText:


    def test_text_shows_both_strand_labels(self, view):
        text = view.to_text(60)
        assert "5'" in text and "3'" in text

    def test_text_contains_the_construct(self, view, plan):
        text = view.to_text(60)
        assert plan.construct_forward[:40] in text.replace("\n", "")

    def test_text_has_no_trailing_whitespace(self, view):
        for line in view.to_text(60).splitlines():
            assert line == line.rstrip()

    def test_narrow_width_still_pairs(self, view):
        rows = view.rows(10)
        assert "".join(row.ticks for row in rows) == view.paired()


class TestFragmentView:


    def test_built_from_the_oligos_not_from_the_construct(self, plan):
        for fragment in plan.fragments:
            single = fragment_duplex(fragment)
            assert single.top.rstrip(GAP) == fragment.forward
            assert single.bottom.lstrip(GAP) == complement(
                reverse_complement(fragment.reverse)
            )
            assert single.mismatches() == []

    def test_each_fragment_shows_its_two_overhangs(self, plan):
        for fragment in plan.fragments:
            single = fragment_duplex(fragment)
            left, right = len(fragment.left_overhang), len(fragment.right_overhang)
            assert single.top[:left] == fragment.left_overhang
            assert single.bottom[:left] == GAP * left
            assert single.top[-right:] == GAP * right
            assert single.bottom[-right:] == complement(fragment.right_overhang)

    def test_fragment_views_line_up_with_the_construct_view(self, plan, view):

        for span, fragment in zip(view.top_fragments, plan.fragments, strict=False):
            single = fragment_duplex(fragment)
            assert view.top[span.start : span.end] == single.top.strip(GAP)

        for span, fragment in zip(view.bottom_fragments, plan.fragments, strict=False):
            single = fragment_duplex(fragment)
            assert view.bottom[span.start : span.end] == single.bottom.strip(GAP)

    def test_the_stagger_is_the_same_in_both_views(self, plan, view):

        for fragment, top_span, bottom_span in zip(
            plan.fragments, view.top_fragments, view.bottom_fragments
        , strict=False):
            assert bottom_span.start - top_span.start == fragment.bottom_offset


class TestJunctionView:


    def build(self, pair: str = "NdeI / XhoI", insert_length: int = 60):
        from gsynth_engine import vectors
        from gsynth_engine.cloning import clone
        from gsynth_engine.tests.test_cloning import build_vector, clean_filler

        left, right = (p.strip() for p in pair.split("/"))
        vector = (
            vectors.sequence_of("pET-21a")["sequence"]
            if pair == "NdeI / XhoI" else build_vector(left, right)
        )
        design = design_extended_sequence(
            clean_filler(insert_length, 3), enzyme_pair=pair,
            target_oligo_length=200,
        )
        return clone(
            vector, design.construct_forward,
            insert_reverse=design.construct_reverse,
            left_enzyme=left, right_enzyme=right, name="X",
        )

    def view_of(self, result, index: int = 0):
        junction = result.junctions[index]
        return junction_view(
            result.plasmid, name=junction.name, enzyme=junction.enzyme,
            position=junction.position, overhang=junction.overhang,
            kind=junction.kind, strand="", flank=14,
        )

    def test_the_joined_strands_pair_everywhere(self):
        view = self.view_of(self.build())
        assert view.joined_pairs.count("|") == len(view.joined_top)

    def test_the_joined_top_is_really_the_plasmid(self):

        result = self.build()
        junction = result.junctions[0]
        view = self.view_of(result)

        expected = "".join(
            result.plasmid[(junction.position + i) % len(result.plasmid)]
            for i in range(-14, 14)
        )
        assert view.joined_top == expected

    def test_both_pieces_carry_the_overhang_on_opposite_strands(self):

        view = self.view_of(self.build())
        low, high = view.overhang_span


        left_top_bases = len(view.left_top.rstrip(GAP))
        left_bottom_bases = len(view.left_bottom.rstrip(GAP))
        assert left_top_bases != left_bottom_bases

        right_top_gap = len(view.right_top) - len(view.right_top.lstrip(GAP))
        right_bottom_gap = len(view.right_bottom) - len(view.right_bottom.lstrip(GAP))
        assert right_top_gap != right_bottom_gap
        assert high - low == len(view.overhang)

    def test_the_overhang_span_holds_the_enzymes_overhang(self):
        view = self.view_of(self.build())
        low, high = view.overhang_span
        assert view.joined_top[low:high] == view.overhang
        assert view.compatible

    def test_a_three_prime_overhang_staggers_the_other_way(self):

        five = self.view_of(self.build("NdeI / XhoI"))
        three = self.view_of(self.build("KpnI / SacI"))

        assert five.kind == "5'" and three.kind == "3'"

        assert five.overhang_span[0] >= five.seam
        assert three.overhang_span[1] <= three.seam

    def test_a_blunt_junction_has_no_overhang_to_show(self):
        view = self.view_of(self.build("EcoRV / SmaI"))
        assert view.kind == "blunt"
        assert view.overhang_span[0] == view.overhang_span[1]
        assert view.compatible

    def test_both_junctions_of_a_construct_are_drawable(self):
        result = self.build()
        for index in range(len(result.junctions)):
            view = self.view_of(result, index)
            assert view.compatible, view.reason
            assert view.width == 28

    def test_a_corrupted_seam_is_reported_not_drawn_as_fine(self):

        result = self.build()
        junction = result.junctions[0]
        at = junction.position
        broken = (
            result.plasmid[:at]
            + ("G" if result.plasmid[at] != "G" else "C")
            + result.plasmid[at + 1:]
        )
        view = junction_view(
            broken, name=junction.name, enzyme=junction.enzyme,
            position=junction.position, overhang=junction.overhang,
            kind=junction.kind, strand="", flank=14,
        )
        assert not view.compatible
        assert junction.overhang in view.reason

    def test_an_empty_plasmid_is_refused(self):
        with pytest.raises(SequenceError, match="plasmid is empty"):
            junction_view("", name="x", enzyme="NdeI", position=0,
                          overhang="TA", kind="5'", strand="")
