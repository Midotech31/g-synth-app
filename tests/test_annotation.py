"""Tests for plasmid feature auto-annotation."""
from gsynth_core.annotation import annotate_features, KNOWN_FEATURES
from gsynth_core.io.records import SeqRecord


def test_detects_t7_promoter():
    seq = "AAA" + "TAATACGACTCACTATAGGG" + "AAA"
    rec = SeqRecord(sequence=seq, name="test")
    hits = annotate_features(rec)
    assert any(h.name == "T7_promoter" for h in hits)


def test_detects_6xhis_tag():
    seq = "ATG" + "CATCATCATCATCATCAT" + "TAA"
    rec = SeqRecord(sequence=seq, name="test")
    hits = annotate_features(rec)
    assert any(h.name == "6xHis_tag" for h in hits)


def test_detects_reverse_strand():
    """Feature on reverse strand should be detected with strand=-1."""
    from gsynth_core.sequence import reverse_complement
    pat = KNOWN_FEATURES["T7_promoter"][1]
    seq = "AAA" + reverse_complement(pat) + "AAA"
    rec = SeqRecord(sequence=seq, name="test")
    hits = annotate_features(rec)
    rev = [h for h in hits if h.name == "T7_promoter" and h.strand == -1]
    assert rev


def test_append_to_record():
    seq = "TAATACGACTCACTATAGGG"
    rec = SeqRecord(sequence=seq, name="test")
    annotate_features(rec, append_to_record=True)
    assert any(f.qualifiers.get("label") == "T7_promoter" for f in rec.features)
