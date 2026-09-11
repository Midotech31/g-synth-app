import io
from dataclasses import replace

import pytest
from Bio import SeqIO

from gsynth_engine import vectors
from gsynth_engine.align import Scoring, align
from gsynth_engine.duplex import DuplexView, Span
from gsynth_engine.gel import (
    GEL_LADDERS,
    ladder_payload,
    recommended_ladder,
    restriction_digest_sizes,
)
from gsynth_engine.genbank import to_genbank
from gsynth_engine.hybridization import hybridize
from gsynth_engine.ligation import nanograms, plan_ligation, ratio_of
from gsynth_engine.preflight import ssd_preflight
from gsynth_engine.provenance import build_provenance
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.thermo import BufferConditions, duplex_thermodynamics, melting_temperature


def test_empty_gc_and_alignment_similarity():
    assert gc_content('') == 0
    assert align('AAAA', 'CCCC', mode='local').similarity == 0
    assert Scoring(matrix={'A': {'A': 5}}).score('A', 'N') == -4


def test_duplex_reports_mismatch_coordinates():
    view = DuplexView(top='ACGT', bottom='TGGC')
    assert view.paired() == '||xx'
    assert view.mismatches() == [2, 3]
    assert Span('region', 2, 7).length == 5


def test_digest_topology_and_ladder_payload():
    dna = 'AAAAGAATTCCCCCCC'
    assert restriction_digest_sizes(dna, ['EcoRI']) == [len(dna)]
    assert restriction_digest_sizes(dna, ['EcoRI'], circular=False) == [11, 5]
    assert restriction_digest_sizes(dna, ['EcoRI', 'EcoRI'], circular=False) == [11, 5]
    assert recommended_ladder([]) == 'broad-range'
    assert {p['key']: tuple(p['bands']) for p in ladder_payload()} == {
        key: entry['bands'] for key, entry in GEL_LADDERS.items()
    }
    with pytest.raises(SequenceError, match='at least one'):
        restriction_digest_sizes(dna, [])


def test_genbank_multiline_comments_round_trip():
    note = 'Source identifier and sequence review. ' * 8
    record = SeqIO.read(io.StringIO(to_genbank('ATG', comments=[note, ''])), 'genbank')
    assert ' '.join(record.annotations['comment'].split()) == note.strip()


@pytest.mark.parametrize('length', [0, -1])
def test_mass_conversion_rejects_invalid_lengths(length):
    with pytest.raises(SequenceError, match='positive'):
        nanograms(1, length)


def test_mass_total_and_undefined_ratio():
    plan = plan_ligation(vector_length=1000, insert_length=100, vector_ng=10, ratio=2)
    assert plan.total_ng == 12
    with pytest.raises(SequenceError, match='greater than zero'):
        ratio_of(vector_length=1000, insert_length=100, vector_ng=0, insert_ng=1)


@pytest.mark.parametrize('coding,sequence', [(False, 'GCT' * 4), (True, 'ATG' + 'GCT' * 80)])
def test_ssd_preflight_records_length_and_start(coding, sequence):
    result = design_small_sequence(sequence, is_coding=coding)
    report = ssd_preflight(result)
    assert report.workflow == 'ssd'
    assert report.verdict in {'ready', 'review'}
    assert {check.code for check in report.checks} >= {'SSD_START_CODON'}


@pytest.mark.parametrize('reads,count', [({}, 0), ({'read': 'ACGT'}, 1), (None, None)])
def test_provenance_read_payload_is_summarized(reads, count):
    result = build_provenance('verify', parameters={'reads': reads}, output_sequence='ACGT')
    assert result['parameters']['reads'] == {'count': count}


def test_monovalent_limit_is_continuous():
    zero = BufferConditions(name='zero', oligo_nM=250, na_mM=100, mg_mM=0)
    trace = replace(zero, mg_mM=0.0001)
    assert melting_temperature('ACGTTGCAAGGCTTAGCCAT', conditions=zero) == pytest.approx(
        melting_temperature('ACGTTGCAAGGCTTAGCCAT', conditions=trace)
    )
    with pytest.raises(KeyError, match='parameters'):
        duplex_thermodynamics('AN')


def test_missing_and_corrupt_vector_assets(monkeypatch, tmp_path, request):
    vectors.sequence_of.cache_clear()
    request.addfinalizer(vectors.sequence_of.cache_clear)
    monkeypatch.setattr(vectors, 'DATA', tmp_path)
    spec = vectors.DEFAULT_VECTOR
    assert vectors.sequence_of(spec.key) is None
    vectors.sequence_of.cache_clear()
    (tmp_path / spec.bundled).write_text('{"sequence":"ATG"}')
    with pytest.raises(ValueError, match='bundled sequence'):
        vectors.sequence_of(spec.key)
    assert replace(spec, tags=()).tag_summary == 'no tags'
    assert 'C-terminal His-tag' in spec.tag_summary


def test_unknown_optional_vector_enzyme_is_ignored():
    record = vectors.sequence_of(vectors.DEFAULT_VECTOR.key)
    spec = replace(vectors.DEFAULT_VECTOR, unique_sites=(*vectors.DEFAULT_VECTOR.unique_sites, 'unknown'))
    assert vectors.validate(record['sequence'], spec).matches


def test_hybridization_display_preserves_both_strands():
    result = hybridize('AATTATGC', 'GGCCGCAT')
    rows = result.rows(3)
    assert ''.join(row['top'] for row in rows) == result.top
    assert ''.join(row['bottom'] for row in rows) == result.bottom
    assert result.paired_percent == 100
    assert rows[-1]['bottom_end'] == 1
    assert replace(result, overlap_end=result.overlap_start).paired_percent == 0


def test_hybridization_states_and_ambiguity():
    result = hybridize('ATGCCGTA', 'TACGGCAT', analysis_temperature_c=-100)
    assert result.predicted_state == 'favourable_at_temperature'
    assert replace(result, analysis_temperature_c=150).predicted_state == 'temperature_above_tm'
    assert replace(result, tm_c=None).predicted_state == 'not_scored'
    assert replace(result, mismatches=1).predicted_state == 'mismatches_not_thermodynamically_scored'
    repeated = hybridize('AAAAAA', 'TTTT')
    assert repeated.alternative_placements == 2
    assert any('ambiguous' in warning for warning in repeated.warnings)
    with pytest.raises(SequenceError, match='too large'):
        hybridize('A' * 2001, 'T' * 2001)
