"""Tests for sequence operations and translation."""
import pytest

from gsynth_core.sequence import (
    clean_dna, complement, find_orfs, gc_content, hamming_distance,
    is_palindrome, reverse_complement, six_frame_translate, translate,
    validate_dna, ensure_dna, iupac_match, find_motif, find_motif_both_strands,
)


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATGC") == "GCAT"
        assert reverse_complement("AAAA") == "TTTT"

    def test_iupac_full(self):
        """G-Synth 2.x bug fix: full IUPAC complement preserved."""
        assert reverse_complement("ATGCRYSWATG") == "CATWSRYGCAT"
        assert reverse_complement("RYSWKMBDHVN") == "NBDHVKMWSRY"

    def test_case_preserved(self):
        assert reverse_complement("atgc") == "gcat"
        assert reverse_complement("AtGc") == "gCaT"

    def test_empty(self):
        assert reverse_complement("") == ""

    def test_palindrome(self):
        assert is_palindrome("GAATTC")  # EcoRI
        assert is_palindrome("GGATCC")  # BamHI
        assert is_palindrome("ATGCAT")  # NsiI — IS palindromic (RC = self)
        assert not is_palindrome("ATGCC")  # asymmetric


class TestGCContent:
    def test_50_pct(self):
        assert gc_content("ATGC") == 50.0
        assert gc_content("AAGG") == 50.0

    def test_100_pct(self):
        assert gc_content("GGGCCC") == 100.0

    def test_0_pct(self):
        assert gc_content("AAATTT") == 0.0

    def test_as_fraction(self):
        assert gc_content("ATGC", as_fraction=True) == 0.5


class TestCleanValidate:
    def test_clean_strips_whitespace(self):
        assert clean_dna("ATG \n CGT") == "ATGCGT"

    def test_clean_strict_drops_iupac(self):
        assert clean_dna("ATGN", allow_ambiguous=False) == "ATG"

    def test_clean_keeps_iupac(self):
        assert clean_dna("ATGN", allow_ambiguous=True) == "ATGN"

    def test_validate_empty(self):
        r = validate_dna("")
        assert not r.ok

    def test_validate_min_length(self):
        r = validate_dna("ATGC", min_length=10)
        assert not r.ok

    def test_ensure_raises(self):
        with pytest.raises(ValueError):
            ensure_dna("")


class TestIUPACMatch:
    def test_n_matches_any(self):
        assert iupac_match("N", "A")
        assert iupac_match("NNN", "ATG")

    def test_r_matches_ag(self):
        assert iupac_match("R", "A")
        assert iupac_match("R", "G")
        assert not iupac_match("R", "C")

    def test_unequal_length(self):
        assert not iupac_match("ATG", "AT")


class TestFindMotif:
    def test_basic(self):
        assert find_motif("AAAATGAATG", "ATG") == [3, 7]

    def test_iupac(self):
        assert find_motif("AAGAATTCAA", "RAATTY") == [2]

    def test_both_strands_palindrome(self):
        hits = find_motif_both_strands("AAAGAATTCAAA", "GAATTC")
        assert hits["+"] == [3]
        assert hits["-"] == []   # palindrome, RC == self


class TestHamming:
    def test_basic(self):
        assert hamming_distance("ATGC", "ATGC") == 0
        assert hamming_distance("ATGC", "ATGT") == 1

    def test_mismatch_length(self):
        with pytest.raises(ValueError):
            hamming_distance("ATG", "ATGC")


class TestTranslate:
    def test_basic(self):
        assert translate("ATGGCT") == "MA"

    def test_stop_codon(self):
        assert translate("ATGTAATGA") == "M**"

    def test_to_stop(self):
        assert translate("ATGGCTTAATGA", to_stop=True) == "MA"

    def test_frame(self):
        # frame=1 skips the first nt → reads ATG GCT from positions 1-6
        assert translate("AATGGCT", frame=1) == "MA"

    def test_six_frame(self):
        result = six_frame_translate("ATGGCTGAA")
        assert "+1" in result and "-1" in result
        assert result["+1"][0] == "M"


class TestFindORFs:
    def test_simple_orf(self):
        orfs = find_orfs("ATGGCTAAA" + "TAA", min_length_aa=2)
        assert len(orfs) >= 1
        assert orfs[0].strand == "+"
        assert orfs[0].protein.startswith("M")

    def test_finds_reverse_strand(self):
        """G-Synth 2.x bug: only forward strand was searched."""
        # Build a seq with an ORF on the reverse strand
        rev_orf = "ATGCGTGCT" + "TAA"
        seq = "AAAAAA" + reverse_complement(rev_orf) + "AAAAAA"
        orfs = find_orfs(seq, min_length_aa=2, both_strands=True)
        rev_orfs = [o for o in orfs if o.strand == "-"]
        assert len(rev_orfs) >= 1

    def test_min_length_filter(self):
        seq = "ATGGCT" + "TAA"  # 2-aa ORF
        assert find_orfs(seq, min_length_aa=10) == []
