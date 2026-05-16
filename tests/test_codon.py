"""Tests for CAI and codon optimization."""
import pytest
from gsynth_core.codon import (
    calculate_cai, codon_frequencies, optimize, OptimizationParams,
    relative_adaptiveness, reverse_translate_best,
)
from gsynth_core.constants import CODON_TABLE_ECOLI
from gsynth_core.sequence import translate


class TestCAI:
    def test_perfect_cai(self):
        """A sequence using all most-frequent codons should give CAI = 1.0."""
        protein = "MKRAILVFDGSYE*"
        best = "".join(max(CODON_TABLE_ECOLI[a], key=CODON_TABLE_ECOLI[a].__getitem__)
                       for a in protein)
        cai = calculate_cai(best, organism="e_coli")
        assert cai > 0.99

    def test_low_cai_worst(self):
        """Worst codons give low CAI."""
        protein = "KRADENF"
        worst = "".join(min(CODON_TABLE_ECOLI[a], key=CODON_TABLE_ECOLI[a].__getitem__)
                        for a in protein)
        cai = calculate_cai(worst, organism="e_coli")
        assert cai < 0.5

    def test_unknown_organism(self):
        with pytest.raises(ValueError):
            calculate_cai("ATG", organism="alien")

    def test_relative_adaptiveness_max_is_1(self):
        w = relative_adaptiveness("e_coli")
        # For every AA, at least one codon should have w = 1.0
        for aa in "ACDEFGHIKLMNPQRSTVWY*":
            codons = [c for c in CODON_TABLE_ECOLI[aa]]
            assert any(abs(w[c] - 1.0) < 1e-6 for c in codons)


class TestOptimize:
    def test_round_trip_preserves_protein(self):
        protein = "MKRLAVFVLLFLAWLG"
        res = optimize(protein, params=OptimizationParams(organism="e_coli"))
        assert translate(res.optimized_dna) == protein

    def test_removes_forbidden_sites(self):
        protein = "MKRLAVFVLLFLAWLGAAQCKRESCGIVHTNSCDK"
        res = optimize(protein, params=OptimizationParams(
            organism="e_coli",
            avoid_sites=("GAATTC", "GGATCC", "AAGCTT"),
        ))
        assert not res.forbidden_sites_remaining

    def test_high_cai(self):
        protein = "MKRLAVFVLLFLAWLG"
        res = optimize(protein, params=OptimizationParams(organism="e_coli"))
        assert res.cai_after > 0.85

    def test_no_homopolymer_above_threshold(self):
        protein = "MKKKKKKKKKKK"  # would normally give all AAA
        res = optimize(protein, params=OptimizationParams(
            organism="e_coli", max_homopolymer=5,
        ))
        import re
        assert not re.search(r"(.)\1{5,}", res.optimized_dna)

    def test_empty_protein_raises(self):
        with pytest.raises(ValueError):
            optimize("")
