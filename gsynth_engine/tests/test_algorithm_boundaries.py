import random
from dataclasses import replace

import pytest

from gsynth_engine import annotations, cloning, codon, esd, pcr, vectors
from gsynth_engine.genbank import Feature
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.tests.test_cloning import build_vector, clean_filler
from gsynth_engine.tests.test_preflight import GENE


@pytest.mark.parametrize('sequence', ['', 'ACGTZ'])
def test_detector_rejects_invalid_input(sequence):
    with pytest.raises(SequenceError, match='IUPAC'):
        annotations.detect_common_features(sequence)


def test_invalid_existing_coordinates_do_not_hide_motif():
    motif = annotations.COMMON_MOTIFS[0]
    found = annotations.detect_common_features(motif.sequence, existing=[{
        'name': motif.name, 'type': motif.feature_type, 'start': None, 'end': 'bad', 'direction': 1,
    }])
    assert any(item['annotation']['name'] == motif.name for item in found)


def test_palindromic_motif_is_reported_once(monkeypatch):
    motif = replace(annotations.COMMON_MOTIFS[0], sequence='ACGT')
    monkeypatch.setattr(annotations, 'COMMON_MOTIFS', (motif,))
    assert len(annotations.detect_common_features('ACGT', circular=False)) == 1


def test_manual_annotation_basis_exports_as_note():
    assert Feature.from_dict({'basis': 'Source record'}).qualifiers['note'] == 'Source record'


@pytest.mark.parametrize('counts,message', [({}, 'missing='), ({key: 0 for key in codon._CODON_TO_AA}, 'No observations')])
def test_invalid_codon_counts_fail(counts, message):
    with pytest.raises(ValueError, match=message):
        codon._weights_from_counts('test', counts)


def test_rare_codon_constraint_and_empty_scoring_window():
    rare = next(iter(codon.ECOLI.rare()))
    constraints = codon.Constraints(avoid_rare=True)
    assert any('low-frequency' in item[2] for item in codon._violations(rare, constraints, codon.ECOLI))
    assert codon._cost('ATG', constraints, codon.ECOLI, window=(3, 3)) == 0
    with pytest.raises(SequenceError, match='empty'):
        codon.optimise('TAA')


def test_perturbation_preserves_protein():
    codons = ['GCT']
    codon._perturb(codons, 0, 3, 'A', codon.ECOLI, random.Random(0))
    assert codons != ['GCT']
    assert codon.translate(''.join(codons)) == 'A'
    codon._perturb(codons, 100, 3, 'A', codon.ECOLI, random.Random(0))
    assert codon.translate(''.join(codons)) == 'A'


def test_fragment_with_one_site_cannot_use_same_enzyme_twice():
    with pytest.raises(SequenceError, match='only once'):
        cloning.digest_linear('AAAAGAATTCCCC', left_enzyme='EcoRI', right_enzyme='EcoRI')


def test_blunt_clone_without_reverse_and_summary():
    result = cloning.clone(build_vector('EcoRV', 'SmaI'), 'ATGGCTTAA', left_enzyme='EcoRV', right_enzyme='SmaI', orf_start=0)
    assert result.gc == round(gc_content(result.plasmid), 1)
    assert result.is_clonable
    assert cloning.open_reading_frames('ATG') == []
    assert cloning._motif_of(vectors.DEFAULT_VECTOR, 'absent') == ''


def test_annotation_with_removed_coding_region_retains_only_flank():
    vector = build_vector('NdeI', 'XhoI')
    backbone = cloning.linearise(vector, left_enzyme='NdeI', right_enzyme='XhoI')
    start = (backbone.vector_start + backbone.length - 6) % len(vector)
    feature = {'name': 'flank', 'type': 'CDS', 'start': start, 'end': start + 15,
        'translation_start': start + 9, 'translation_end': start + 15}
    moved = cloning._remap_annotations([feature], backbone, backbone.length + 20)[0]
    assert moved['type'] == 'misc_feature'
    assert 'translation_start' not in moved and 'translation_end' not in moved
    assert cloning._remap_annotations([feature], replace(backbone, circular_source=False), 10) == []


@pytest.mark.parametrize('spacing,status', [(2, 'review'), (8, 'pass'), (30, 'block')])
def test_frame_context_uses_regulatory_class_and_distance(spacing, status):
    report = cloning._assess_reading_frame('C' * 60 + 'ATGGCTTAA' + 'C' * 31,
        insert_start=60, insert_end=69, orf_start=0, protein='MA', stop_at=66,
        annotations=[{'name': 'initiation', 'type': 'regulatory', 'regulatory_class': 'ribosome_binding_site',
            'start': 54 - spacing, 'end': 60 - spacing, 'direction': 1}],
        vector_spec=None, tags=[], existing_problems=[], start_source='declared')
    assert next(check for check in report.checks if check.code == 'FRAME_RBS_CONTEXT').status == status


def test_no_stop_frame_is_blocked():
    result = cloning.clone('CATATG' + 'CCC' * 10 + 'CTCGAG' + 'CCC' * 20,
        'TATG' + 'CCC' * 8, left_enzyme='NdeI', right_enzyme='XhoI', orf_start=1)
    assert any('No stop codon' in warning for warning in result.warnings)
    assert any(check.status == 'block' for check in result.reading_frame.checks if check.code == 'FRAME_TRANSLATED_TERMINUS')


def test_small_primer_utilities_and_display():
    assert pcr._clamp_of(0) == ''
    assert pcr._self_complementarity('ATG') == 0
    assert pcr._cross_dimer('ATG', 'CAT') == 0
    result = pcr.design_pcr(GENE)
    assert result.insert == GENE
    assert result.forward.as_row['Length (nt)'] == len(result.forward.sequence)
    with pytest.raises(SequenceError, match='not enough template'):
        pcr._pick_anneal('ATG', start=0, direction=1)


@pytest.mark.parametrize('primer,message', [('', 'both primer'), ('A' * 301, '300'), ('CCC' + GENE[:24], 'unpaired')])
def test_invalid_custom_primers_are_rejected(primer, message):
    with pytest.raises(SequenceError, match=message):
        pcr._resolve_custom_primer(primer, template=GENE, boundary=0, direction=1, allow_tail=False)


def test_edited_reverse_primer_requires_selected_site():
    result = pcr.design_pcr(GENE, left_enzyme='NdeI', right_enzyme='XhoI')
    with pytest.raises(SequenceError, match="reverse primer's"):
        pcr.design_pcr(GENE, left_enzyme='NdeI', right_enzyme='XhoI',
            forward_primer=result.forward.sequence, reverse_primer=result.reverse.anneals)


def test_short_restriction_flanks_are_reported():
    result = pcr.design_pcr(GENE, left_enzyme='NdeI', right_enzyme='XhoI')
    edited = pcr.design_pcr(GENE, left_enzyme='NdeI', right_enzyme='XhoI',
        forward_primer=result.forward.sequence[4:], reverse_primer=result.reverse.sequence[4:])
    assert any('fewer than six' in warning for warning in edited.warnings)
    assert edited.insert == edited.digest.top


def test_assembly_terminal_collision_is_reported():
    pool = esd._OverhangPool({'AAGC'})
    assert 'terminal restriction' in pool.problem('AAGC')


def test_assembly_fragment_count_adapts_to_overhang_size():
    plan = esd.design_extended_sequence(clean_filler(100), target_oligo_length=24, overhang_length=8)
    assert plan.verify() == []


def test_ssd_start_validation_precedes_stop_removal():
    with pytest.raises(SequenceError, match='begin with ATG'):
        design_small_sequence('TAA', is_coding=True, remove_stop=True)
