"""Tests for CRISPR scoring and guide design."""
import pytest
from gsynth_core.crispr import (
    cfd_score, mit_specificity_score, doench_2016_score,
    design_guides, CasType, find_pam_sites,
)


class TestCFD:
    def test_perfect_match(self):
        g = "GTCACCTCCAATGACTAGGG"
        assert cfd_score(g, g, pam="TGG") > 0.99

    def test_bad_pam(self):
        g = "GTCACCTCCAATGACTAGGG"
        # TAA is not in the penalty table → 0
        assert cfd_score(g, g, pam="TAA") == 0.0

    def test_mismatch_near_pam_heavy_penalty(self):
        g = "GTCACCTCCAATGACTAGGG"
        # Mismatch at position 20 (very close to PAM)
        off = g[:19] + "A"
        score_near = cfd_score(g, off, pam="AGG")
        assert score_near < 0.1

    def test_mismatch_far_from_pam_light_penalty(self):
        g = "GTCACCTCCAATGACTAGGG"
        # Mismatch at position 1
        off = "A" + g[1:]
        score_far = cfd_score(g, off, pam="AGG")
        assert score_far > 0.5

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            cfd_score("ATG", "ATGCATGC", pam="AGG")

    def test_wrong_length_raises(self):
        with pytest.raises(ValueError):
            cfd_score("AAA", "AAA")


class TestDoench:
    def test_TTTT_kills_score(self):
        """TTTT is a Pol III terminator → activity drops."""
        # 30-mer = 4 upstream + 20 protospacer + 3 PAM + 3 downstream
        # Context WITH TTTT in protospacer
        ctx_with_tttt = "AAAA" + "ATTTTCCAATGACTAGGGT" + "T" + "TGG" + "TTT"
        assert len(ctx_with_tttt) == 30, f"len = {len(ctx_with_tttt)}"
        ctx_without = "AAAA" + "GTCACCTCCAATGACTAGGG" + "TGG" + "TTT"
        assert len(ctx_without) == 30, f"len = {len(ctx_without)}"
        assert doench_2016_score(ctx_with_tttt) < doench_2016_score(ctx_without)

    def test_score_in_unit_range(self):
        ctx = "AAAA" + "GTCACCTCCAATGACTAGGG" + "TGG" + "TTT"
        s = doench_2016_score(ctx)
        assert 0.0 <= s <= 1.0


class TestPAMSearch:
    def test_finds_ngg(self):
        seq = "CCAAATGGTACGGCTTGG"
        hits = find_pam_sites(seq, CasType.SpCas9)
        # NGG occurs at TGG (pos 5), CGG (pos 11), TGG (pos 16)
        assert any(p == 5 and s == "+" for p, s in hits)


class TestDesignGuides:
    def test_finds_guides_in_gene(self):
        gene = ("ATGAAACGCCTGGCTGTTTTTGTGCTGCTGTTCCTGGCTTGGCTGGGCGCTGCTCAGCGCAGCGCA"
                "GAATTCGCGCGCCGCTGACTAGGGCCAGGCCAAGGTCACCTCCAATGACTAGGGTGGTTATG")
        guides = design_guides(gene, cas=CasType.SpCas9, min_on_target=0.3)
        assert len(guides) > 0
        for g in guides:
            assert 0 <= g.on_target <= 1.0
            assert g.cas == CasType.SpCas9

    def test_off_target_annotation(self):
        gene = ("ATGAAACGCCTGGCTGTTTTTGTGCTGCTGTTCCTGGCTTGGCTGGGCGCTGCTCAGCGCAGCGCA"
                "GAATTCGCGCGCCGCTGACTAGGGCCAGGCCAAGGTCACCTCCAATGACTAGGGTGGTTATG")
        guides = design_guides(gene, cas=CasType.SpCas9, min_on_target=0.3,
                               genome_reference=gene * 3, top_n=3)
        for g in guides:
            assert g.specificity is not None
            assert g.off_target_count is not None
