"""Tests for I/O modules."""
import pytest
from gsynth_core.io import (
    parse_fasta, write_fasta, parse_genbank, write_genbank,
    SeqRecord, Feature, FeatureLocation,
)


class TestFASTA:
    def test_parse_single(self):
        recs = parse_fasta(">seq1\nATGCATGC\n")
        assert len(recs) == 1
        assert recs[0].name == "seq1"
        assert recs[0].sequence == "ATGCATGC"

    def test_parse_multi(self):
        recs = parse_fasta(">a desc1\nATGC\n>b\nGGGGCCCC\n")
        assert len(recs) == 2
        assert recs[0].description == "desc1"

    def test_wrapped(self):
        recs = parse_fasta(">long\nATGCATGC\nATGCATGC\n")
        assert recs[0].sequence == "ATGCATGCATGCATGC"

    def test_round_trip(self):
        original = ">a\nATGCATGC\n>b\nGGGGCCCC\n"
        recs = parse_fasta(original)
        out = write_fasta(recs)
        recs2 = parse_fasta(out)
        assert recs[0].sequence == recs2[0].sequence
        assert recs[1].sequence == recs2[1].sequence


class TestGenBank:
    def test_parse_and_write(self):
        gb = """LOCUS       test                     100 bp    DNA     circular SYN 01-JAN-2025
DEFINITION  test plasmid
FEATURES             Location/Qualifiers
     CDS             10..60
                     /label="orf"
                     /codon_start=1
ORIGIN
        1 atgaaacgcc tggctgtttt tgtgctgctg ttcctggctt ggttgggcgc tgctcagcgc
       61 agcgcagaat tcgcgcgccg ctgactaggg ccaggccaag
//
"""
        recs = parse_genbank(gb)
        assert len(recs) == 1
        assert recs[0].topology == "circular"
        assert len(recs[0].features) == 1
        assert recs[0].features[0].type == "CDS"

        # Round-trip
        out = write_genbank(recs[0])
        recs2 = parse_genbank(out)
        assert recs[0].sequence.upper() == recs2[0].sequence.upper()


class TestRecordTypes:
    def test_seqrecord_basic(self):
        r = SeqRecord(sequence="ATGC", name="test")
        assert len(r) == 4

    def test_add_feature(self):
        r = SeqRecord(sequence="ATGCATGCATGC", name="test")
        f = Feature(type="CDS",
                    location=FeatureLocation(start=0, end=3, strand=1),
                    qualifiers={"label": "start"})
        r.add_feature(f)
        assert len(r.features) == 1
        assert f.name == "start"

    def test_features_in(self):
        r = SeqRecord(sequence="A" * 100, name="test")
        r.add_feature(Feature("CDS", FeatureLocation(10, 30)))
        r.add_feature(Feature("CDS", FeatureLocation(50, 60)))
        assert len(r.features_in(0, 40)) == 1
        assert len(r.features_in(40, 70)) == 1
        assert len(r.features_in(70, 100)) == 0
