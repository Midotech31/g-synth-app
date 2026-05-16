"""Tests for pairwise alignment and MSA."""
import pytest
from gsynth_core.alignment import align_global, align_local, multiple_alignment


class TestPairwise:
    def test_identical(self):
        aln = align_global("ATGCATGC", "ATGCATGC")
        assert aln.identity_percent == 100.0

    def test_classic_rosalind(self):
        """GATTACA vs GCATGCU is a classic textbook example."""
        aln = align_global("GATTACA", "GCATGCU")
        assert aln.score > 0
        assert "GATTACA" in aln.aligned_a.replace("-", "")

    def test_local_finds_motif(self):
        aln = align_local("AAAAAGATTACAAAAA", "TTTTTGATTACATTTT")
        assert "GATTACA" in aln.aligned_a.replace("-", "")
        assert aln.identity_percent == 100.0

    def test_local_no_match(self):
        aln = align_local("AAAA", "GGGG")
        assert aln.score == 0.0


class TestMSA:
    def test_star_alignment(self):
        seqs = [("a", "ATGAAACGCCTGGCTGTTTTT"),
                ("b", "ATGAAACGTCTGGCAGTTTTT"),
                ("c", "ATGAAACGTCTCGCTGTTTTT")]
        msa = multiple_alignment(seqs, method="star")
        # All aligned sequences must have same length
        L = len(msa.aligned[0])
        assert all(len(s) == L for s in msa.aligned)

    def test_singleton(self):
        msa = multiple_alignment([("only", "ATGC")])
        assert len(msa.aligned) == 1

    def test_empty(self):
        msa = multiple_alignment([])
        assert len(msa.aligned) == 0
