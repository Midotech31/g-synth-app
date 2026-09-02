"""Tests for antiparallel strand hybridization and cohesive-end geometry."""
from __future__ import annotations

import pytest

from gsynth_engine.hybridization import hybridize
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.thermo import BufferConditions


class TestAntiparallelPlacement:
    def test_two_ordered_strands_form_a_flush_duplex(self):
        first = "ATGCCGTA"
        result = hybridize(first, reverse_complement(first))

        assert result.complementarity == "exact"
        assert result.mismatches == 0
        assert result.paired_bases == len(first)
        assert result.top == first
        assert result.bottom == "TACGGCAT"
        assert set(result.marks) == {"|"}
        assert result.left_end["kind"] == "blunt"
        assert result.right_end["kind"] == "blunt"

    def test_the_two_external_five_prime_sticky_ends_are_visible(self):
        # Physical duplex:
        # 5′ AATTATGC     3′
        #        TACGCCGG 5′
        result = hybridize("AATTATGC", "GGCCGCAT")

        assert result.top == "AATTATGC    "
        assert result.bottom == "    TACGCCGG"
        assert result.marks == "    ||||    "
        assert result.overlap_length == 4
        assert result.left_end == {
            "end": "left", "strand": "first", "polarity": "5′",
            "sequence": "AATT", "length": 4, "start": 0, "end_position": 4,
        }
        assert result.right_end["strand"] == "second"
        assert result.right_end["polarity"] == "5′"
        assert result.right_end["sequence"] == "GGCC"

    def test_three_prime_extensions_are_not_mislabelled_as_five_prime(self):
        result = hybridize("ATGCAATT", "GCATGGCC")

        assert result.left_end["strand"] == "second"
        assert result.left_end["polarity"] == "3′"
        assert result.left_end["sequence"] == "GGCC"
        assert result.right_end["strand"] == "first"
        assert result.right_end["polarity"] == "3′"
        assert result.right_end["sequence"] == "AATT"

    def test_an_internal_mismatch_remains_visible_and_disables_tm(self):
        first = "ATGCCGTA"
        partner = list(reverse_complement(first))
        partner[3] = "A" if partner[3] != "A" else "C"
        result = hybridize(first, "".join(partner))

        assert result.mismatches == 1
        assert "×" in result.marks
        assert result.complementarity == "partial"
        assert result.tm_c is None
        assert any("mismatched" in warning for warning in result.warnings)


class TestThermodynamicScope:
    def test_tm_uses_the_conditions_named_in_the_result(self):
        conditions = BufferConditions(
            name="custom anneal", oligo_nM=20_000, na_mM=75, mg_mM=1.5,
        )
        result = hybridize(
            "ATGCCGTAGCTAGCTA",
            reverse_complement("ATGCCGTAGCTAGCTA"),
            conditions=conditions,
            analysis_temperature_c=37,
        )

        assert result.tm_c is not None
        assert result.tm_margin_c == round(result.tm_c - 37, 1)
        assert result.conditions.summary == "20 µM total strand, 75 mM Na⁺, 1.5 mM Mg²⁺"
        assert result.delta_h_kcal_mol is not None
        assert result.delta_s_cal_mol_k is not None

    def test_less_than_four_pairs_is_not_called_a_cohesive_end(self):
        result = hybridize("AAA", "TTT")
        assert result.complementarity == "insufficient"
        assert result.predicted_state == "insufficient_complementarity"
        assert any("fewer than 4" in warning for warning in result.warnings)


class TestInput:
    def test_invalid_dna_is_refused(self):
        with pytest.raises(SequenceError, match="not A, C, G or T"):
            hybridize("ATUG", "CAT")

    def test_empty_strands_are_refused(self):
        with pytest.raises(SequenceError, match="first strand is empty"):
            hybridize("", "ATGC")

    def test_unreasonable_temperature_is_refused(self):
        with pytest.raises(SequenceError, match="temperature"):
            hybridize("ATGC", "GCAT", analysis_temperature_c=200)
