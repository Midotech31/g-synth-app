"""Tests for codon optimisation.

One property carries the module: **the protein never changes**. A codon
optimiser that improves expression and quietly alters a residue has done
something far worse than nothing, and the mistake is invisible until the
protein misbehaves. It is asserted here for every sequence the suite
touches, including randomly generated ones.
"""
from __future__ import annotations

import random

import pytest

from gsynth_engine.cloning import find_sites, translate
from gsynth_engine.codon import (
    CODON_DATA_SHA256,
    CODON_DATA_VERSION,
    ECOLI,
    SYNONYMS,
    TABLES,
    Constraints,
    back_translate,
    build_table,
    codon_adaptation_index,
    optimise,
)
from gsynth_engine.sequence import SequenceError, gc_content

# An enterocin-like gene: AT-rich, full of codons E. coli reads slowly, and
# carrying an internal NdeI site that would break NdeI/XhoI cloning.
DONOR_GENE = (
    "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA"
    "TTAAATTTAGCATTATTAGGATTAGCAAGTTTATTAGGTAAAGGTATTAGTAAATTAGGA"
)
DONOR_PROTEIN = "MTTSKLGKGLGYIGNNGAHMGLNLALLGLASLLGKGISKLG"


def random_protein(length: int, seed: int) -> str:
    rng = random.Random(seed)
    residues = [aa for aa in SYNONYMS if aa != "*"]
    return "M" + "".join(rng.choice(residues) for _ in range(length - 1))


# ── The invariant ───────────────────────────────────────────────────────────


class TestProteinIsPreserved:
    def test_the_donor_gene_keeps_its_protein(self):
        result = optimise(DONOR_GENE)
        assert result.protein == DONOR_PROTEIN
        assert translate(result.sequence).rstrip("*") == DONOR_PROTEIN

    def test_holds_for_random_proteins(self):
        """Including the awkward ones: single-codon residues, rare families."""
        for seed in range(12):
            protein = random_protein(60, seed)
            result = optimise(protein, is_protein=True)
            assert translate(result.sequence).rstrip("*") == protein, seed

    def test_holds_under_every_constraint_at_once(self):
        constraints = Constraints(
            avoid_enzymes=("NdeI", "XhoI", "BamHI", "EcoRI"),
            avoid_motifs=("AAAAAA", "GGGGGG"),
            max_homopolymer=4,
            gc_min=40.0,
            gc_max=60.0,
        )
        for seed in range(6):
            protein = random_protein(80, seed)
            result = optimise(protein, is_protein=True, constraints=constraints)
            assert translate(result.sequence).rstrip("*") == protein, seed

    def test_length_follows_from_the_protein(self):
        result = optimise(DONOR_PROTEIN, is_protein=True, keep_stop=False)
        assert result.length == 3 * len(DONOR_PROTEIN)

        with_stop = optimise(DONOR_PROTEIN, is_protein=True, keep_stop=True)
        assert with_stop.length == 3 * len(DONOR_PROTEIN) + 3
        assert translate(with_stop.sequence).endswith("*")

    def test_auto_preserves_a_non_methionine_mature_peptide(self):
        peptide = "GIVEQCCTSICSLYQLENYCG"
        result = optimise(peptide, is_protein=True, keep_stop=False)
        assert result.input_protein == peptide
        assert result.protein == peptide
        assert result.protein_context == "mature_peptide"
        assert result.initiator_methionine_added is False
        assert result.recommended_design_is_coding is False
        assert translate(result.sequence) == peptide

    def test_auto_detects_an_existing_initiator_methionine(self):
        result = optimise("MGIVEQ", is_protein=True, keep_stop=False)
        assert result.protein == "MGIVEQ"
        assert result.protein_context == "complete_orf"
        assert result.initiator_methionine_added is False
        assert result.recommended_design_is_coding is True
        assert result.sequence.startswith("ATG")

    def test_complete_orf_adds_exactly_one_initiator_methionine(self):
        result = optimise(
            "GIVEQ", is_protein=True, protein_context="complete_orf", keep_stop=False,
        )
        assert result.input_protein == "GIVEQ"
        assert result.protein == "MGIVEQ"
        assert result.initiator_methionine_added is True
        assert result.recommended_design_is_coding is True
        assert translate(result.sequence) == "MGIVEQ"

    def test_mature_override_preserves_a_peptide_that_starts_with_methionine(self):
        result = optimise(
            "MGIVEQ", is_protein=True, protein_context="mature_peptide", keep_stop=False,
        )
        assert result.protein == "MGIVEQ"
        assert result.protein_context == "mature_peptide"
        assert result.recommended_design_is_coding is False


# ── What optimisation is for ────────────────────────────────────────────────


class TestAdaptation:
    def test_rare_codons_are_removed(self):
        low_frequency_gene = "ATG" + "AGG" * 12 + "CTA" * 12
        result = optimise(low_frequency_gene)
        assert result.rare_codons_before == 24
        assert result.rare_codons_after == 0

    def test_the_adaptation_index_rises(self):
        result = optimise(DONOR_GENE)
        assert result.cai_before is not None
        assert result.cai_after > result.cai_before
        assert result.cai_after > 0.8

    def test_gc_moves_into_the_window(self):
        """The donor is AT-rich; synthesis and PCR both want it nearer 50%."""
        result = optimise(DONOR_GENE)
        assert result.gc_before is not None and result.gc_before < 35
        assert 40 <= result.gc_after <= 60

    def test_it_reports_how_much_it_changed(self):
        result = optimise(DONOR_GENE)
        assert 0 < result.changed_codons <= len(DONOR_PROTEIN)


# ── Constraints ─────────────────────────────────────────────────────────────


class TestConstraints:
    def test_the_cloning_sites_are_removed(self):
        """A gene with an internal NdeI site cannot be cloned NdeI/XhoI,
        however well it translates."""
        assert find_sites(DONOR_GENE, "NdeI", circular=False)

        result = optimise(
            DONOR_GENE, constraints=Constraints(avoid_enzymes=("NdeI", "XhoI"))
        )
        assert find_sites(result.sequence, "NdeI", circular=False) == []
        assert find_sites(result.sequence, "XhoI", circular=False) == []
        assert "CATATG" in result.sites_removed
        assert result.warnings == []

    def test_arbitrary_motifs_are_removed(self):
        result = optimise(
            DONOR_PROTEIN, is_protein=True,
            constraints=Constraints(avoid_motifs=("GGTGGT", "CTGCTG")),
        )
        assert "GGTGGT" not in result.sequence
        assert "CTGCTG" not in result.sequence

    def test_homopolymers_are_kept_short(self):
        from gsynth_engine.sequence import longest_homopolymer

        result = optimise(
            random_protein(120, 3), is_protein=True,
            constraints=Constraints(max_homopolymer=4),
        )
        assert longest_homopolymer(result.sequence) <= 4

    def test_local_gc_is_checked_in_a_window(self):
        """A sequence can average 50% and still have an unsynthesisable patch."""
        result = optimise(
            random_protein(200, 5), is_protein=True,
            constraints=Constraints(gc_min=40, gc_max=60, gc_window=50),
        )
        sequence = result.sequence
        for start in range(0, len(sequence) - 50 + 1, 3):
            local = gc_content(sequence[start : start + 50])
            assert 40 <= local <= 60, f"{local:.0f}% at {start}"

    def test_unsatisfiable_constraints_are_reported_not_hidden(self):
        """Trp and Met have one codon each, so the GC cannot be moved at all."""
        result = optimise(
            "MWWWWWWWWWWMM", is_protein=True,
            constraints=Constraints(gc_min=10, gc_max=20),
        )
        assert result.warnings, "a GC window it cannot reach must be stated"

    def test_a_site_that_cannot_be_removed_blocks_the_design(self):
        """Poly-tryptophan has one codon, so TGGTGG cannot be avoided."""
        result = optimise(
            "MWWWWWWWW", is_protein=True,
            constraints=Constraints(avoid_motifs=("TGGTGG",)),
        )
        assert not result.is_clean
        assert any("TGGTGG" in p for p in result.problems)

    def test_severity_follows_the_consequence(self):
        """A leftover rare codon costs translation speed; a leftover site
        costs the whole cloning strategy. Calling both problems would train
        the user to ignore the word."""
        result = optimise(
            DONOR_GENE,
            constraints=Constraints(
                avoid_enzymes=("NdeI", "XhoI"), gc_min=45, gc_max=55,
            ),
        )
        # The site is gone, so nothing blocks the design...
        assert result.is_clean, result.problems
        assert find_sites(result.sequence, "NdeI", circular=False) == []
        # ...even if a compromise elsewhere is worth mentioning.
        assert isinstance(result.warnings, list)

    def test_a_remaining_enzyme_site_is_not_duplicated_as_a_warning(self):
        result = optimise(
            "SM",
            is_protein=True,
            constraints=Constraints(avoid_enzymes=("CviAII",)),
            max_rounds=0,
        )
        assert any("CATG" in problem for problem in result.problems)
        assert result.warnings == []


# ── Determinism ─────────────────────────────────────────────────────────────


class TestDeterminism:
    def test_the_same_input_gives_the_same_gene(self):
        """Someone has ordered this sequence; it must be reproducible."""
        first = optimise(DONOR_GENE, constraints=Constraints(avoid_enzymes=("NdeI",)))
        second = optimise(DONOR_GENE, constraints=Constraints(avoid_enzymes=("NdeI",)))
        assert first.sequence == second.sequence

    def test_a_different_seed_may_differ_but_still_encodes_the_protein(self):
        constraints = Constraints(gc_min=45, gc_max=55)
        a = optimise(random_protein(100, 9), is_protein=True, constraints=constraints, seed=1)
        b = optimise(random_protein(100, 9), is_protein=True, constraints=constraints, seed=2)
        assert translate(a.sequence).rstrip("*") == translate(b.sequence).rstrip("*")


# ── The usage table ─────────────────────────────────────────────────────────


class TestCodonTable:
    def test_every_codon_has_a_weight(self):
        assert len(ECOLI.weights) == 64

    def test_each_amino_acid_has_one_codon_at_full_weight(self):
        """That is what relative adaptiveness means."""
        for amino_acid, codons in SYNONYMS.items():
            top = max(ECOLI.weight(c) for c in codons)
            assert top == pytest.approx(1.0), amino_acid

    def test_every_bundled_host_is_complete_and_normalized(self):
        assert set(TABLES) == {
            "ecoli", "b_subtilis", "p_putida", "l_lactis",
            "c_glutamicum", "s_coelicolor", "s_cerevisiae", "k_phaffii",
            "k_lactis", "y_lipolytica", "h_sapiens", "c_griseus",
            "s_frugiperda", "d_melanogaster", "n_benthamiana",
        }
        for key, table in TABLES.items():
            assert len(table.weights) == 64, key
            assert table.taxon_id is not None, key
            assert table.dataset in {"RefSeq", "GenBank"}, key
            assert table.dataset_release == "September 2021", key
            assert table.coding_sequences and table.coding_sequences > 0, key
            assert table.codon_count and table.codon_count > 0, key
            assert table.gc_percent and 0 < table.gc_percent < 100, key
            for amino_acid, codons in SYNONYMS.items():
                assert max(table.weight(codon) for codon in codons) == pytest.approx(1.0), (
                    key, amino_acid
                )

    def test_host_selection_changes_codons_without_changing_protein(self):
        protein = "MLRLKFY"
        bacterial = back_translate(protein, table=TABLES["ecoli"])
        yeast = back_translate(protein, table=TABLES["s_cerevisiae"])
        assert bacterial != yeast
        assert translate(bacterial) == protein
        assert translate(yeast) == protein

    def test_host_profiles_are_numerically_distinct_not_labels(self):
        codons = sorted(ECOLI.weights)
        signatures = {
            tuple(round(table.weights[codon], 12) for codon in codons)
            for table in TABLES.values()
        }
        assert len(signatures) == len(TABLES)

    @pytest.mark.parametrize(("host", "amino_acid", "preferred"), [
        ("ecoli", "R", "CGC"),
        ("b_subtilis", "K", "AAA"),
        ("p_putida", "K", "AAG"),
        ("l_lactis", "A", "GCT"),
        ("c_glutamicum", "K", "AAG"),
        ("s_coelicolor", "A", "GCC"),
        ("s_cerevisiae", "L", "TTG"),
        ("k_phaffii", "L", "TTG"),
        ("k_lactis", "F", "TTC"),
        ("y_lipolytica", "A", "GCC"),
        ("h_sapiens", "F", "TTC"),
        ("c_griseus", "I", "ATC"),
        ("s_frugiperda", "Y", "TAC"),
        ("d_melanogaster", "I", "ATC"),
        ("n_benthamiana", "A", "GCT"),
    ])
    def test_species_profiles_preserve_published_rankings(
        self, host: str, amino_acid: str, preferred: str,
    ):
        assert TABLES[host].best(amino_acid) == preferred

    def test_the_known_e_coli_preferences(self):
        """Pinned preferences are independently visible in the source rows."""
        assert ECOLI.best("L") == "CTG"
        assert ECOLI.best("R") == "CGC"
        assert ECOLI.best("I") == "ATT"
        assert ECOLI.best("*") == "TAA"

    def test_the_notorious_rare_codons_are_rare(self):
        """The 0.1 threshold is applied to the committed species profile."""
        rare = ECOLI.rare()
        assert rare == {"AGG", "CTA"}

    def test_the_source_snapshot_is_pinned(self):
        assert CODON_DATA_VERSION == "September 2021"
        assert CODON_DATA_SHA256 == (
            "4c6d35b4c42dd449b9e441e6eef9cb17777211036f0cbc8e7c5f27616b9a8ab0"
        )

    def test_single_codon_residues_are_never_rare(self):
        rare = ECOLI.rare()
        assert "ATG" not in rare and "TGG" not in rare

    def test_a_table_can_be_built_from_reference_genes(self):
        """The honest way to get a CAI: measure it against genes you chose."""
        table = build_table(["ATGCTGCTGCTGAAA", "ATGCTGCTGTTAAAA"], name="test")
        assert table.best("L") == "CTG"
        assert table.weight("TTA") == pytest.approx(1 / 5)
        assert table.name == "test"

    def test_an_empty_reference_set_is_refused(self):
        with pytest.raises(SequenceError, match="No complete codons"):
            build_table([""], name="empty")


class TestAdaptationIndex:
    def test_all_preferred_codons_scores_one(self):
        best = back_translate("MLRIKAG")
        assert codon_adaptation_index(best) == pytest.approx(1.0)

    def test_rare_codons_score_low(self):
        assert codon_adaptation_index("CTAAGGCTAAGG") < 0.1

    def test_single_codon_residues_are_excluded(self):
        """Met and Trp offer no choice, so counting them inflates every gene."""
        assert codon_adaptation_index("ATGTGG") == 0.0

    def test_it_is_measured_against_the_table_given(self):
        table = build_table(["TTATTATTA"], name="odd")
        assert codon_adaptation_index("TTATTA", table) == pytest.approx(1.0)
        assert codon_adaptation_index("CTGCTG", table) < 1.0


# ── Input handling ──────────────────────────────────────────────────────────


class TestInput:
    def test_a_partial_codon_is_refused_with_the_count(self):
        with pytest.raises(SequenceError, match="not a multiple of three"):
            optimise("ATGAAAA")

    def test_invalid_bases_are_refused(self):
        with pytest.raises(SequenceError, match="not A, C, G or T"):
            optimise("ATGXYZATG")

    def test_invalid_residues_are_refused(self):
        with pytest.raises(SequenceError, match="not amino acids"):
            optimise("MKBJZ", is_protein=True)

    def test_an_empty_protein_is_refused(self):
        with pytest.raises(SequenceError):
            optimise("", is_protein=True)

    def test_an_unknown_protein_context_is_refused(self):
        with pytest.raises(SequenceError, match="Protein context"):
            optimise("MGIVEQ", is_protein=True, protein_context="unknown")

    def test_a_trailing_stop_is_not_translated_twice(self):
        with_stop = optimise(DONOR_GENE + "TAA", keep_stop=True)
        assert with_stop.protein == DONOR_PROTEIN
        assert translate(with_stop.sequence).count("*") == 1

    def test_back_translate_takes_the_preferred_codon(self):
        assert back_translate("ML") == "ATGCTG"

    def test_back_translate_refuses_a_bad_residue(self):
        with pytest.raises(SequenceError, match="not amino acids"):
            back_translate("MKZ")
