from __future__ import annotations

import math

import pytest

from gsynth_engine.sequence import melting_temperature as sequence_tm
from gsynth_engine.thermo import (
    ANNEALING,
    PRIMER,
    BufferConditions,
    duplex_thermodynamics,
    melting_temperature,
)


def test_enthalpy_and_entropy_match_the_published_table_for_a_poly_a():

    enthalpy, entropy = duplex_thermodynamics("AAAAAAAA")

    expected_h = 7 * -7.9 + 2 * 2.3
    expected_s = 7 * -22.2 + 2 * 4.1

    assert enthalpy == pytest.approx(expected_h)
    assert entropy == pytest.approx(expected_s)


def test_palindrome_carries_the_symmetry_penalty():

    enthalpy, entropy = duplex_thermodynamics("GCGCGCGC")

    expected_h = 4 * -9.8 + 3 * -10.6 + 2 * 0.1
    expected_s = 4 * -24.4 + 3 * -27.2 + 2 * -2.8 - 1.4

    assert enthalpy == pytest.approx(expected_h)
    assert entropy == pytest.approx(expected_s)


def test_a_dimer_and_its_complement_share_parameters():

    assert duplex_thermodynamics("CA") == pytest.approx(duplex_thermodynamics("TG"))
    assert duplex_thermodynamics("GT") == pytest.approx(duplex_thermodynamics("AC"))


def test_duplex_thermodynamics_is_stable_and_ordered():

    enthalpy, entropy = duplex_thermodynamics("GTAAAACGACGGCCAGT")

    assert enthalpy < 0
    assert entropy < 0


def test_tm_depends_on_base_order_not_only_on_gc_content():

    sequences = ["GGCCGGAATTCC", "GAGCGAGCTAGC", "CCGGCCGGAATT"]
    assert len({len(s) for s in sequences}) == 1
    assert len({sum(1 for b in s if b in "GC") for s in sequences}) == 1

    temperatures = [melting_temperature(s) for s in sequences]
    assert len({round(t, 1) for t in temperatures}) == len(sequences)


def test_gc_rich_melts_higher_than_at_rich_at_equal_length():
    assert melting_temperature("GCGCGCGCGC") > melting_temperature("ATATATATAT")


def test_tm_rises_with_length():
    short = melting_temperature("GTAAAACGAC")
    long = melting_temperature("GTAAAACGACGGCCAGT")
    assert long > short


def test_one_molar_sodium_is_the_uncorrected_reference_state():

    sequence = "GTAAAACGACGGCCAGT"
    enthalpy, entropy = duplex_thermodynamics(sequence)
    raw_kelvin = (enthalpy * 1000.0) / (entropy + 1.9872 * math.log(500e-9 / 4))

    assert melting_temperature(sequence, na_mM=1000.0) == pytest.approx(
        raw_kelvin - 273.15, abs=1e-6
    )


def test_lower_salt_lowers_tm():
    sequence = "GTAAAACGACGGCCAGT"
    assert melting_temperature(sequence, na_mM=1000.0) > melting_temperature(
        sequence, na_mM=50.0
    )


def test_magnesium_raises_tm():

    sequence = "GTAAAACGACGGCCAGT"
    without = melting_temperature(sequence, na_mM=50.0)
    with_mg = melting_temperature(sequence, na_mM=50.0, mg_mM=1.5)

    assert with_mg > without + 3.0


def test_dntps_chelate_magnesium():

    sequence = "GTAAAACGACGGCCAGT"
    chelated = melting_temperature(sequence, na_mM=50.0, mg_mM=1.5, dntp_mM=1.5)
    no_magnesium = melting_temperature(sequence, na_mM=50.0)

    assert chelated == pytest.approx(no_magnesium)


def test_m13_primer_lands_in_the_published_range():

    tm = melting_temperature("GTAAAACGACGGCCAGT", na_mM=50.0, mg_mM=1.5)
    assert 55.0 <= tm <= 60.0


def test_concentration_shifts_tm():
    sequence = "GTAAAACGACGGCCAGT"
    assert melting_temperature(sequence, oligo_nM=1000.0) > melting_temperature(
        sequence, oligo_nM=100.0
    )


def test_self_complementary_strand_uses_the_undivided_concentration_term():

    sequence = "GGCCAATTGGCC"
    enthalpy, entropy = duplex_thermodynamics(sequence)
    expected = (enthalpy * 1000.0) / (
        entropy + 1.9872 * math.log(500e-9)
    ) - 273.15

    assert melting_temperature(sequence, na_mM=1000.0) == pytest.approx(
        expected, abs=1e-6
    )


def test_short_sequences_fall_back_to_wallace():

    assert melting_temperature("GATC") == pytest.approx(12.0)
    assert melting_temperature("GGCC") == pytest.approx(16.0)
    assert melting_temperature("AAAA") == pytest.approx(8.0)


def test_eight_nucleotides_is_the_first_nearest_neighbour_length():

    wallace = melting_temperature("GCGCGCG")
    assert wallace == pytest.approx(2 * 0 + 4 * 7)

    thermodynamic = melting_temperature("GCGCGCGC")
    assert thermodynamic != pytest.approx(4 * 8)


def test_empty_sequence_is_zero():
    assert melting_temperature("") == 0.0
    assert duplex_thermodynamics("") == (0.0, 0.0)


def test_annealing_conditions_match_the_protocol_reaction():

    assert ANNEALING.oligo_nM == 50_000.0
    assert ANNEALING.na_mM == 50.0
    assert ANNEALING.mg_mM == 0.0


def test_the_annealing_reaction_melts_higher_than_a_primer_dilution():

    sequence = "GTAAAACGACGGCCAGT"
    assert melting_temperature(sequence, conditions=ANNEALING) > melting_temperature(
        sequence, conditions=PRIMER
    ) + 5.0


def test_keywords_override_individual_values_within_conditions():
    sequence = "GTAAAACGACGGCCAGT"
    assert melting_temperature(sequence, conditions=ANNEALING, mg_mM=1.5) > (
        melting_temperature(sequence, conditions=ANNEALING)
    )
    assert melting_temperature(
        sequence, conditions=ANNEALING, oligo_nM=PRIMER.oligo_nM
    ) == pytest.approx(melting_temperature(sequence, conditions=PRIMER))


def test_conditions_summarise_themselves_for_a_report():
    plain = BufferConditions(name="test", oligo_nM=500.0, na_mM=50.0)
    assert plain.summary == "0.5 µM total strand, 50 mM Na⁺"

    with_salts = BufferConditions(
        name="pcr", oligo_nM=500.0, na_mM=50.0, mg_mM=1.5, dntp_mM=0.8
    )
    assert "1.5 mM Mg²⁺" in with_salts.summary
    assert "0.8 mM dNTP" in with_salts.summary


def test_sequence_module_delegates_to_the_nearest_neighbour_model():

    sequence = "GTAAAACGACGGCCAGT"
    assert sequence_tm(sequence) == melting_temperature(sequence)
    assert sequence_tm(sequence, mg_mM=1.5) == melting_temperature(
        sequence, mg_mM=1.5
    )


def test_input_is_cleaned_before_calculation():

    plain = melting_temperature("GTAAAACGACGGCCAGT")
    messy = melting_temperature(">primer\ngtaaaacg acggcc\nAGT\n")
    assert messy == pytest.approx(plain)
