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
        item for item in detect_common_features("CCCAAGGAGAAAAATGCCC", existing=[
            {"name": "target", "type": "CDS", "start": 13, "end": 19, "direction": 1},
        ])
        if item["annotation"]["type"] == "RBS"
    )

    assert match["annotation"]["name"] == "Shine-Dalgarno candidate"
    assert match["matched_sequence"] == "AAGGAG"
    assert "review" in match["basis"].lower()


def test_sd_requires_a_downstream_possible_start_in_transcript_direction():
    sequence = "CCCAAGGAG" + "C" * 7 + "ATG" + "C" * 20
    forward = [m for m in detect_common_features(sequence) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]
    reverse = [m for m in detect_common_features(reverse_complement(sequence)) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]
    assert len(forward) == len(reverse) == 1
    assert forward[0]['annotation']['direction'] == 1
    assert reverse[0]['annotation']['direction'] == -1
    assert '7 nt upstream' in forward[0]['basis']
    assert 'mRNA' in forward[0]['basis']
    assert forward[0]['annotation']['inferred'] is True
    for isolated in ['AAGGAG' + 'C' * 30, 'ATG' + 'C' * 7 + 'AAGGAG', 'AAGGAGATG']:
        assert not [m for m in detect_common_features(isolated) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]


def test_sd_overlapping_spellings_do_not_make_two_annotations():
    matches = detect_common_features('AAGGAGG' + 'C' * 7 + 'ATG')
    assert len([m for m in matches if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]) == 1


def test_sd_circular_context_and_coding_overlap_are_respected():
    sequence = 'AAGGAG' + 'C' * 7 + 'ATGCCC'
    rotated = sequence[10:] + sequence[:10]
    assert len([m for m in detect_common_features(rotated, circular=True) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]) == 1
    assert not [m for m in detect_common_features(sequence, existing=[
        {'name': 'coding region', 'type': 'CDS', 'start': 0, 'end': len(sequence), 'direction': 1},
    ]) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]


def test_existing_rbs_alias_and_strand_do_not_create_a_false_duplicate():
    sequence = 'AAGGAG' + 'C' * 7 + 'ATG'
    known = {'name': 'RBS', 'type': 'RBS', 'start': 0, 'end': 6, 'direction': 1}
    assert not [m for m in detect_common_features(sequence, existing=[known]) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]
    known['direction'] = -1
    assert len([m for m in detect_common_features(sequence, existing=[known]) if ('Shine-Dalgarno' in m['annotation']['name'] or 'SD-like' in m['annotation']['name'])]) == 1


def test_t7_terminator_direction_follows_transcription_in_the_reference():
    from gsynth_engine.vectors import sequence_of
    sequence = sequence_of('pET-21a')['sequence']
    terminators = [m for m in detect_common_features(sequence, circular=True) if m['annotation']['type'] == 'terminator']
    assert len(terminators) == 1
    assert terminators[0]['annotation']['direction'] == -1


def test_unassigned_sd_spelling_is_not_exported_as_a_regulatory_rbs():
    match = detect_common_features('AAGGAG' + 'C' * 7 + 'ATG')[0]
    assert match['annotation']['name'] == 'SD-like motif (unassigned)'
    assert match['annotation']['type'] == 'misc_feature'
    assert 'regulatory_class' not in match['annotation']
    assert 'No annotated CDS' in match['basis']
    assert 'coordinate 14' in match['basis']


def test_sd_prefers_the_annotated_translation_start_over_an_incidental_codon():
    sequence = 'AAGGAGCCCCTTGCCATGCCCCCC'
    match = detect_common_features(sequence, existing=[
        {'name': 'target', 'type': 'CDS', 'start': 13, 'end': 23, 'direction': 1,
         'translation_start': 15, 'translation_end': 23},
    ])[0]
    assert match['annotation']['type'] == 'RBS'
    assert '9 nt upstream' in match['basis']
    assert 'annotated CDS: target' in match['basis']
    assert 'coordinate 16' in match['basis']


def test_reverse_and_circular_sd_link_to_the_same_annotated_start():
    sequence = 'AAGGAGCCCCCCCATGCCC'
    reverse = reverse_complement(sequence)
    match = detect_common_features(reverse, existing=[
        {'name': 'reverse target', 'type': 'CDS', 'start': 0, 'end': 6, 'direction': -1},
    ])[0]
    assert match['annotation']['type'] == 'RBS'
    assert 'reverse target' in match['basis']
    assert 'coordinate 6 on the reverse strand' in match['basis']
    rotated = sequence[10:] + sequence[:10]
    match = detect_common_features(rotated, circular=True, existing=[
        {'name': 'wrapped target', 'type': 'CDS', 'start': 3, 'end': 9, 'direction': 1},
    ])[0]
    assert match['annotation']['type'] == 'RBS'
    assert 'wrapped target' in match['basis']


def test_overlapping_sd_spellings_prefer_the_cds_associated_match():
    sequence = 'AAGGAGG' + 'CCCCTTG' + 'C' * 7 + 'ATGCCC'
    matches = detect_common_features(sequence, existing=[
        {'name': 'target', 'type': 'CDS', 'start': 21, 'end': len(sequence), 'direction': 1},
    ])
    sd = [m for m in matches if 'SD-like' in m['annotation']['name'] or m['annotation']['type'] == 'RBS']
    assert len(sd) == 1
    assert sd[0]['annotation']['type'] == 'RBS'
    assert sd[0]['annotation']['start'] == 1
    assert '14 nt upstream' in sd[0]['basis']


def test_pet_recombinant_keeps_one_correctly_oriented_terminator_and_the_original_rbs():
    from gsynth_engine.cloning import clone
    from gsynth_engine.ssd import design_small_sequence
    from gsynth_engine.vectors import sequence_of

    for key in ['pET-21', 'pET-21a']:
        record = sequence_of(key)
        known = next(a for a in record['annotations'] if a['name'] == 'T7 terminator')
        assert known['direction'] == -1
        assert not [m for m in detect_common_features(record['sequence'], circular=True,
                    existing=record['annotations']) if m['annotation']['type'] == 'terminator']
    record = sequence_of('pET-21a')
    ssd = design_small_sequence('GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA',
                                enzyme_pair='NdeI / XhoI')
    result = clone(record['sequence'], ssd.forward, insert_reverse=ssd.reverse,
                   left_enzyme='NdeI', right_enzyme='XhoI', orf_start=ssd.orf_start,
                   vector_annotations=record['annotations'])
    assert len(result.plasmid) == 5490
    assert result.reading_frame.rbs_start == 5352
    assert result.reading_frame.rbs_spacing_nt == 8
    assert next(a for a in result.annotations if a['name'] == 'T7 terminator')['direction'] == 1
    matches = detect_common_features(result.plasmid, circular=True, existing=result.annotations)
    assert not [m for m in matches if m['annotation']['type'] in ('terminator', 'RBS')]
    orphan = next(m for m in matches if m['annotation']['start'] == 141)
    assert orphan['annotation']['type'] == 'misc_feature'
    assert 'unassigned' in orphan['annotation']['name']
