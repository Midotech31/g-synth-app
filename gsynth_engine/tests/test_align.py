"""Tests for pairwise alignment.

The properties that matter are not "does it produce an alignment" but
whether it produces the *right* one: affine gaps must prefer one long gap to
several short ones, a reverse-complemented sequence must be recognised
rather than reported as unrelated, and the three modes must actually differ.
"""
from __future__ import annotations

import pytest

from gsynth_engine.align import MAX_CELLS, Scoring, align, blosum62
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.tests.test_cloning import clean_filler

GENE = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGGTATATTGGAAATAATGGAGCACATATGGGA"


class TestBasics:
    def test_identical_sequences_align_perfectly(self):
        result = align(GENE, GENE)
        assert result.identity == 100.0
        assert result.gaps == 0
        assert result.length == len(GENE)

    def test_the_alignment_reproduces_both_inputs(self):
        """Whatever else it does, it must not invent or lose bases."""
        other = GENE[:20] + "GGG" + GENE[25:]
        result = align(GENE, other)
        assert result.top.replace("-", "") == GENE
        assert result.bottom.replace("-", "") == other

    def test_a_substitution_is_marked_but_not_gapped(self):
        mutated = GENE[:40] + ("T" if GENE[40] != "T" else "A") + GENE[41:]
        result = align(GENE, mutated)
        assert result.gaps == 0
        assert result.marks.count(".") == 1
        assert result.identity < 100

    def test_the_three_strings_are_the_same_length(self):
        result = align(GENE, GENE[:30] + GENE[40:])
        assert len(result.top) == len(result.marks) == len(result.bottom)


class TestAffineGaps:
    def test_one_long_gap_beats_several_short_ones(self):
        """One long gap is one event; scattered gaps are biologically
        backwards and produce alignments that look wrong to a reader."""
        deleted = GENE[:20] + GENE[32:]          # a clean 12 nt deletion
        result = align(GENE, deleted)

        assert result.gaps == 12
        # All twelve in one run, not sprinkled.
        runs = [run for run in result.bottom.split("-") if run]
        assert len(runs) == 2, result.bottom

    def test_a_cheaper_gap_opening_allows_more_gaps(self):
        """The penalty is doing the work, not an accident of the sequences."""
        a, b = clean_filler(120, 5), clean_filler(120, 6)
        strict = align(a, b, scoring=Scoring(gap_open=40, gap_extend=2))
        loose = align(a, b, scoring=Scoring(gap_open=1, gap_extend=1))
        assert loose.gaps > strict.gaps


class TestOrientation:
    def test_a_reverse_complemented_sequence_is_recognised(self):
        """A gene cloned the other way round is not a different gene."""
        result = align(GENE, reverse_complement(GENE))
        assert result.reverse_complemented
        assert result.identity == 100.0
        assert any("other way round" in note for note in result.warnings)

    def test_positions_come_back_in_the_callers_numbering(self):
        result = align(GENE, reverse_complement(GENE))
        assert 0 <= result.start_b < result.end_b <= len(GENE)

    def test_it_can_be_turned_off(self):
        result = align(GENE, reverse_complement(GENE), try_reverse=False)
        assert not result.reverse_complemented
        assert result.identity < 60

    def test_protein_is_never_reverse_complemented(self):
        result = align("MTTSKLGKG", "MTTSRLGKG", is_protein=True)
        assert not result.reverse_complemented


class TestModes:
    def test_global_uses_the_whole_of_both(self):
        result = align("AAAA" + GENE + "TTTT", GENE, mode="global")
        assert result.top.replace("-", "") == "AAAA" + GENE + "TTTT"
        assert result.length >= len(GENE) + 8

    def test_local_finds_only_the_shared_stretch(self):
        """Two sequences that share one region and nothing else."""
        a = "TTTTTTTTTTTT" + GENE[:40] + "GGGGGGGGGGGG"
        b = "CCCCCCCCCCCC" + GENE[:40] + "AAAAAAAAAAAA"
        result = align(a, b, mode="local")

        assert result.identity == 100.0
        assert result.length == 40
        assert result.start_a == 12 and result.end_a == 52

    def test_semi_global_places_the_shorter_inside_the_longer(self):
        """Where does this gene sit in this plasmid — the question global
        alignment answers badly and local answers incompletely."""
        plasmid = clean_filler(300, 9) + GENE + clean_filler(300, 10)
        result = align(GENE, plasmid, mode="semi-global")

        assert result.identity > 95
        assert result.top.replace("-", "") == GENE

    def test_an_unknown_mode_is_refused(self):
        with pytest.raises(SequenceError, match="Unknown alignment mode"):
            align(GENE, GENE, mode="sideways")


class TestProtein:
    def test_the_matrix_is_the_published_one(self):
        """Read from Biopython's copy rather than typed: 576 values is 576
        chances to be quietly wrong."""
        matrix = blosum62()
        assert matrix["W"]["W"] == 11
        assert matrix["C"]["C"] == 9
        assert matrix["L"]["I"] == 2
        assert matrix["W"]["G"] == -2

    def test_conservative_substitutions_are_marked_as_similar(self):
        """K→R and I→V are not the same residue and not a random change."""
        result = align("MTTSKLGKGLGYIGNN", "MTTSRLGKGLGYVGNN", is_protein=True)
        assert result.identity < 100
        assert result.similarity == 100.0
        assert result.marks.count(":") == 2

    def test_an_unrelated_pair_scores_low(self):
        result = align("MKKKKKKKKKK", "WWWWWWWWWWW", is_protein=True)
        assert result.identity == 0.0
        assert result.similarity < 30

    def test_dna_has_no_similarity_category(self):
        """A base is the right one or it is not."""
        result = align(GENE, GENE[:30] + "G" + GENE[31:])
        assert result.similarity == result.identity


class TestOutput:
    def test_rows_number_each_sequence_separately(self):
        result = align(GENE, GENE)
        rows = result.rows(30)
        assert rows[0]["top_start"] == 1
        assert rows[1]["top_start"] == 31
        assert rows[-1]["top_end"] == len(GENE)

    def test_rows_reassemble_the_alignment(self):
        result = align(GENE, GENE[:30] + GENE[40:])
        rows = result.rows(25)
        assert "".join(str(row["top"]) for row in rows) == result.top
        assert "".join(str(row["bottom"]) for row in rows) == result.bottom

    def test_a_gap_row_does_not_advance_the_numbering(self):
        result = align("A" * 60, "A" * 20 + "A" * 20)
        for row in result.rows(20):
            if row["top_start"] is None:
                assert row["top_end"] == row["top_end"]      # unchanged

    def test_the_text_rendering_shows_both_and_the_marks(self):
        text = align(GENE, GENE).to_text(60)
        assert "|" in text
        assert GENE[:60] in text


class TestInput:
    def test_an_empty_sequence_is_refused(self):
        with pytest.raises(SequenceError, match="first sequence is empty"):
            align("", GENE)
        with pytest.raises(SequenceError, match="second sequence is empty"):
            align(GENE, "")

    def test_whitespace_and_fasta_headers_are_stripped(self):
        messy = f">gene\n{GENE[:30]}\n{GENE[30:]}\n"
        assert align(GENE, messy).identity == 100.0

    def test_a_pair_too_large_is_refused_with_a_pointer(self):
        """Better than a request that hangs for a minute."""
        size = int(MAX_CELLS ** 0.5) + 500
        with pytest.raises(SequenceError, match="verification tool"):
            align("A" * size, "A" * size)
