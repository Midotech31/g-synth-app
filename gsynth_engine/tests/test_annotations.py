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

    assert match["annotation"]["name"] == "Shine-Dalgarno candidate"
    assert match["matched_sequence"] == "AAGGAG"
    assert "review" in match["basis"].lower()


def test_sd_requires_a_downstream_possible_start_in_transcript_direction():
    sequence = "CCCAAGGAG" + "C" * 7 + "ATG" + "C" * 20
    forward = [m for m in detect_common_features(sequence) if m['annotation']['type'] == 'RBS']
    reverse = [m for m in detect_common_features(reverse_complement(sequence)) if m['annotation']['type'] == 'RBS']
    assert len(forward) == len(reverse) == 1
    assert forward[0]['annotation']['direction'] == 1
    assert reverse[0]['annotation']['direction'] == -1
    assert '7 nt upstream' in forward[0]['basis']
    assert 'mRNA' in forward[0]['basis']
    assert forward[0]['annotation']['inferred'] is True
    for isolated in ['AAGGAG' + 'C' * 30, 'ATG' + 'C' * 7 + 'AAGGAG', 'AAGGAGATG']:
        assert not [m for m in detect_common_features(isolated) if m['annotation']['type'] == 'RBS']


def test_sd_overlapping_spellings_do_not_make_two_annotations():
    matches = detect_common_features('AAGGAGG' + 'C' * 7 + 'ATG')
    assert len([m for m in matches if m['annotation']['type'] == 'RBS']) == 1


def test_sd_circular_context_and_coding_overlap_are_respected():
    sequence = 'AAGGAG' + 'C' * 7 + 'ATGCCC'
    rotated = sequence[10:] + sequence[:10]
    assert len([m for m in detect_common_features(rotated, circular=True) if m['annotation']['type'] == 'RBS']) == 1
    assert not [m for m in detect_common_features(sequence, existing=[
        {'name': 'coding region', 'type': 'CDS', 'start': 0, 'end': len(sequence), 'direction': 1},
    ]) if m['annotation']['type'] == 'RBS']


def test_existing_rbs_alias_and_strand_do_not_create_a_false_duplicate():
    sequence = 'AAGGAG' + 'C' * 7 + 'ATG'
    known = {'name': 'RBS', 'type': 'RBS', 'start': 0, 'end': 6, 'direction': 1}
    assert not [m for m in detect_common_features(sequence, existing=[known]) if m['annotation']['type'] == 'RBS']
    known['direction'] = -1
    assert len([m for m in detect_common_features(sequence, existing=[known]) if m['annotation']['type'] == 'RBS']) == 1


def test_t7_terminator_direction_follows_transcription_in_the_reference():
    from gsynth_engine.vectors import sequence_of
    sequence = sequence_of('pET-21a')['sequence']
    terminators = [m for m in detect_common_features(sequence, circular=True) if m['annotation']['type'] == 'terminator']
    assert len(terminators) == 1
    assert terminators[0]['annotation']['direction'] == -1
