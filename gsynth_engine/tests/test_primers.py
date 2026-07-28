"""Tests for sequencing primer design.

Run against a real recombinant plasmid rather than a synthetic template,
because the two properties that matter — uniqueness across the whole
molecule, and reading *into* the insert rather than at it — only mean
anything on a molecule with a real backbone in it.
"""
from __future__ import annotations

import pytest

from gsynth_engine import vectors
from gsynth_engine.cloning import clone
from gsynth_engine.primers import (
    DEAD_ZONE,
    SEQUENCING,
    design_sequencing_primers,
)
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.tests.test_cloning import clean_filler
from gsynth_engine.thermo import melting_temperature


def build(insert_length: int, seed: int = 3):
    record = vectors.sequence_of("pET-21a")
    design = design_small_sequence(
        clean_filler(insert_length, seed), enzyme_pair="NdeI / XhoI"
    )
    return clone(
        record["sequence"], design.forward, insert_reverse=design.reverse,
        left_enzyme="NdeI", right_enzyme="XhoI", name="pGS",
    )


@pytest.fixture(scope="module")
def short():
    return build(150)


@pytest.fixture(scope="module")
def long():
    return build(900, seed=4)


class TestFlankingPrimers:
    def test_it_designs_a_primer_on_each_strand(self, short):
        result = design_sequencing_primers(
            result_plasmid := short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
        )
        assert len(result.primers) >= 2
        assert {p.direction for p in result.primers} == {1, -1}
        assert result_plasmid

    def test_every_primer_is_unique_in_the_plasmid(self, short):
        """A primer that binds twice gives a superimposed trace and no data."""
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
        )
        for primer in result.primers:
            forward = short.plasmid.count(primer.sequence)
            reverse = short.plasmid.count(reverse_complement(primer.sequence))
            assert forward + reverse == 1, primer.name

    def test_every_primer_is_really_in_the_template(self, short):
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
        )
        length = len(short.plasmid)
        for primer in result.primers:
            window = "".join(
                short.plasmid[(primer.start + i) % length]
                for i in range(primer.length)
            )
            expected = window if primer.direction == 1 else reverse_complement(window)
            assert primer.sequence == expected, primer.name

    def test_primers_sit_back_from_the_target(self, short):
        """A primer starting at the insert reads its beginning as noise —
        the part with the ATG and the tag in it."""
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
            margin=80,
        )
        forward = [p for p in result.primers if p.direction == 1][0]
        assert forward.reads_from <= short.insert_start

    def test_the_margin_cannot_be_shorter_than_the_dead_zone(self, short):
        with pytest.raises(SequenceError, match="unreadable"):
            design_sequencing_primers(
                short.plasmid,
                target_start=short.insert_start, target_end=short.insert_end,
                margin=DEAD_ZONE - 1,
            )

    def test_melting_temperatures_are_in_the_requested_window(self, short):
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
            tm_min=52, tm_max=62,
        )
        for primer in result.primers:
            assert 52 <= primer.tm <= 62, primer.name

    def test_the_tm_is_the_engines_own_model_under_stated_conditions(self, short):
        """Two numbers in one project must mean the same thing."""
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
        )
        primer = result.primers[0]
        assert primer.tm == pytest.approx(
            round(melting_temperature(primer.sequence, conditions=SEQUENCING), 1)
        )

    def test_lengths_stay_in_range(self, short):
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
            length_range=(20, 24),
        )
        for primer in result.primers:
            assert 20 <= primer.length <= 24


class TestCoverage:
    def test_a_short_insert_is_covered_by_the_two_flanking_primers(self, short):
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
        )
        assert result.covers_target, result.gaps
        assert result.warnings == []

    def test_a_long_insert_gets_internal_primers(self, long):
        """One Sanger read gives ~700 usable bases; a 950 bp insert needs more."""
        result = design_sequencing_primers(
            long.plasmid,
            target_start=long.insert_start, target_end=long.insert_end,
        )
        assert len(result.primers) > 2
        assert result.covers_target, result.gaps

    def test_read_ranges_wrap_the_origin(self, short):
        """A read runs round the origin as readily as any other stretch.
        Clamping instead reported a primer past the origin as reading nothing."""
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
            circular=True,
        )
        assert any(p.reads_to < p.reads_from for p in result.primers), (
            "the insert sits at the end of this plasmid, so a read must wrap"
        )

    def test_an_uncoverable_target_says_so_rather_than_looking_finished(self, long):
        result = design_sequencing_primers(
            long.plasmid,
            target_start=long.insert_start, target_end=long.insert_end,
            read_length=200,
        )
        assert not result.covers_target
        assert any("not reached by any primer" in note for note in result.warnings)

    def test_a_repetitive_insert_reports_what_it_could_not_place(self, short):
        """Nothing in a tandem repeat is unique, so no primer can be made."""
        record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(
            "ATG" + "GCTAGCAAAGGTTTCCGTGAAGATCTGGCAAAATTCCTGCAGGCTAACGGT" * 18,
            enzyme_pair="NdeI / XhoI",
        )
        repetitive = clone(
            record["sequence"], design.forward, insert_reverse=design.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI", name="repeat",
        )
        result = design_sequencing_primers(
            repetitive.plasmid,
            target_start=repetitive.insert_start, target_end=repetitive.insert_end,
        )
        assert any("unique" in note for note in result.warnings)


class TestInput:
    def test_a_target_outside_the_template_is_refused(self, short):
        with pytest.raises(SequenceError, match="does not lie inside"):
            design_sequencing_primers(
                short.plasmid, target_start=0, target_end=len(short.plasmid) + 10,
            )

    def test_an_inverted_target_is_refused(self, short):
        with pytest.raises(SequenceError, match="does not lie inside"):
            design_sequencing_primers(short.plasmid, target_start=500, target_end=100)

    def test_an_empty_template_is_refused(self):
        with pytest.raises(SequenceError, match="template is empty"):
            design_sequencing_primers("", target_start=0, target_end=1)

    def test_the_order_sheet_row_has_what_a_supplier_needs(self, short):
        result = design_sequencing_primers(
            short.plasmid,
            target_start=short.insert_start, target_end=short.insert_end,
        )
        row = result.as_rows[0]
        assert set(row) >= {"Name", "Sequence (5'->3')", "Length (nt)", "Tm (°C)"}
