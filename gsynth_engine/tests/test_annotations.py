"""Common-feature recognition is exact, strand-aware and reviewable."""

from gsynth_engine.annotations import detect_common_features
from gsynth_engine.constants import HIS_TAG
from gsynth_engine.sequence import reverse_complement


def test_detects_known_motifs_on_both_strands():
    sequence = "AAA" + HIS_TAG + "CCC" + reverse_complement(HIS_TAG) + "TTT"
    matches = detect_common_features(sequence)
    his = [match for match in matches if match["annotation"]["name"] == "6×His tag"]

    assert [(item["annotation"]["start"], item["annotation"]["direction"]) for item in his] == [
        (3, 1),
        (24, -1),
    ]
    assert all("Exact match" in item["basis"] for item in his)


def test_detects_a_motif_across_a_circular_origin():
    motif = "TAATACGACTCACTATAGGG"
    sequence = motif[-7:] + "A" * 30 + motif[:-7]
    match = next(
        item for item in detect_common_features(sequence, circular=True)
        if item["annotation"]["name"] == "T7 promoter"
        and item["annotation"]["direction"] == 1
    )

    assert match["annotation"]["start"] == 37
    assert match["annotation"]["end"] == 57
    assert match["annotation"]["end"] > len(sequence)


def test_does_not_duplicate_an_imported_feature_that_covers_the_core():
    motif = "TAATACGACTCACTATAGGG"
    existing = [{"name": "T7 promoter", "start": 0, "end": len(motif) + 2}]

    assert not [
        item for item in detect_common_features(motif + "CC", existing=existing)
        if item["annotation"]["name"] == "T7 promoter"
    ]


def test_does_not_duplicate_a_same_named_import_that_differs_by_one_flanking_base():
    motif = "TAATACGACTCACTATAGGG"
    existing = [{"name": "T7 promoter", "start": 0, "end": len(motif) - 1}]

    assert not [
        item for item in detect_common_features(motif, existing=existing)
        if item["annotation"]["name"] == "T7 promoter"
    ]


def test_linear_sequence_does_not_match_across_its_ends():
    motif = "TAATACGACTCACTATAGGG"
    sequence = motif[-7:] + "A" * 30 + motif[:-7]

    assert not [
        item for item in detect_common_features(sequence)
        if item["annotation"]["name"] == "T7 promoter"
    ]


def test_detects_exact_bacterial_rbs_as_reviewable_evidence():
    match = next(
        item for item in detect_common_features("CCCAAGGAGAAAAATGCCC")
        if item["annotation"]["type"] == "RBS"
    )

    assert match["annotation"]["name"] == "Shine-Dalgarno RBS"
    assert match["matched_sequence"] == "AAGGAG"
    assert "review" in match["basis"].lower()
