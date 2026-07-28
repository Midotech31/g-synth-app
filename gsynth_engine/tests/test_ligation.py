"""Tests for ligation arithmetic.

Short, and constantly got wrong: the quantity that matters is molar and the
quantity you pipette is mass. The tests pin the conversion in both
directions and the advice that surrounds it.
"""
from __future__ import annotations

import pytest

from gsynth_engine.ligation import (
    DALTONS_PER_BP,
    femtomoles,
    ligation_series,
    nanograms,
    plan_ligation,
    ratio_of,
)
from gsynth_engine.sequence import SequenceError


class TestConversion:
    def test_femtomoles_follows_the_published_constant(self):
        """50 ng of a 5,443 bp plasmid, by hand."""
        expected = (50 / (5443 * DALTONS_PER_BP)) * 1e6
        assert femtomoles(50, 5443) == pytest.approx(expected)

    def test_the_two_directions_are_inverses(self):
        for length in (150, 1000, 5443):
            assert nanograms(femtomoles(50, length), length) == pytest.approx(50)

    def test_a_zero_length_is_refused(self):
        with pytest.raises(SequenceError, match="positive number of base pairs"):
            femtomoles(50, 0)


class TestPlan:
    def test_equal_mass_is_nowhere_near_equal_moles(self):
        """The mistake this module exists to prevent."""
        assert ratio_of(
            vector_length=5443, insert_length=150, vector_ng=50, insert_ng=50,
        ) == pytest.approx(36.29, abs=0.05)

    def test_a_three_to_one_ratio_gives_the_expected_mass(self):
        plan = plan_ligation(vector_length=5443, insert_length=150, vector_ng=50, ratio=3)
        # mass_insert = 50 × (150 / 5443) × 3
        assert plan.insert_ng == pytest.approx(50 * 150 / 5443 * 3, abs=0.02)
        assert plan.insert_fmol == pytest.approx(plan.vector_fmol * 3, rel=1e-3)

    def test_the_ratio_comes_back_out(self):
        plan = plan_ligation(vector_length=5443, insert_length=980, vector_ng=50, ratio=5)
        assert ratio_of(
            vector_length=5443, insert_length=980,
            vector_ng=plan.vector_ng, insert_ng=plan.insert_ng,
        ) == pytest.approx(5, abs=0.02)

    def test_sticky_and_blunt_ends_get_different_defaults(self):
        """Blunt ends compete with the vector closing on itself."""
        sticky = plan_ligation(vector_length=5000, insert_length=500, ends="5'")
        blunt = plan_ligation(vector_length=5000, insert_length=500, ends="blunt")
        assert sticky.ratio == 3.0
        assert blunt.ratio == 5.0
        assert blunt.insert_ng > sticky.insert_ng

    def test_blunt_ends_carry_the_advice_that_goes_with_them(self):
        plan = plan_ligation(vector_length=5000, insert_length=500, ends="blunt")
        assert any("dephosphorylate" in note for note in plan.warnings)

    def test_an_unpipettable_amount_is_flagged(self):
        """0.1 ng is below what anyone measures accurately."""
        plan = plan_ligation(vector_length=50_000, insert_length=60, vector_ng=10, ratio=1)
        assert plan.insert_ng < 1
        assert any("pipette accurately" in note for note in plan.warnings)

    def test_too_much_dna_is_flagged(self):
        plan = plan_ligation(vector_length=3000, insert_length=6000, vector_ng=400, ratio=3)
        assert any("a lot for one ligation" in note for note in plan.warnings)

    def test_lengths_the_wrong_way_round_are_noticed(self):
        plan = plan_ligation(vector_length=500, insert_length=5000)
        assert any("right way round" in note for note in plan.warnings)

    def test_it_reports_molecules_when_the_numbers_stop_feeling_real(self):
        plan = plan_ligation(vector_length=5443, insert_length=150, vector_ng=50, ratio=3)
        assert plan.insert_molecules > 1e9

    def test_the_table_names_both_components(self):
        rows = plan_ligation(vector_length=5443, insert_length=150).as_rows()
        assert [row["Component"] for row in rows] == ["vector", "insert"]
        assert "fmol" in rows[0]["Moles"]

    def test_impossible_reactions_are_refused(self):
        with pytest.raises(SequenceError, match="positive number of base pairs"):
            plan_ligation(vector_length=0, insert_length=100)
        with pytest.raises(SequenceError, match="greater than zero"):
            plan_ligation(vector_length=5000, insert_length=100, vector_ng=0)
        with pytest.raises(SequenceError, match="greater than zero"):
            plan_ligation(vector_length=5000, insert_length=100, ratio=0)


class TestSeries:
    def test_it_gives_one_reaction_per_ratio(self):
        """Nobody runs one ligation; they run a small series and pick a plate."""
        series = ligation_series(vector_length=5443, insert_length=150)
        assert [plan.ratio for plan in series] == [1.0, 3.0, 5.0]

    def test_more_insert_at_a_higher_ratio(self):
        series = ligation_series(vector_length=5443, insert_length=150)
        masses = [plan.insert_ng for plan in series]
        assert masses == sorted(masses)

    def test_the_vector_amount_is_the_same_across_the_series(self):
        series = ligation_series(vector_length=5443, insert_length=150, vector_ng=75)
        assert {plan.vector_ng for plan in series} == {75.0}
