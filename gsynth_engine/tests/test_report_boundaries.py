from dataclasses import replace

import pytest

from gsynth_engine import esd, protocol
from gsynth_engine.duplex import DuplexView
from gsynth_engine.ligation import ligation_series
from gsynth_engine.preflight import cloning_preflight
from gsynth_engine.primers import design_sequencing_primers
from gsynth_engine.tests.test_cloning import clean_filler
from gsynth_engine.tests.test_primers import build
from gsynth_engine.thermo import ANNEALING, melting_temperature


@pytest.mark.parametrize('gaps', [[], [(10, 20)]])
def test_worksheet_matches_derived_records(gaps):
    result = build(150)
    primers = design_sequencing_primers(result.plasmid, target_start=result.insert_start, target_end=result.insert_end)
    primers = replace(primers, gaps=gaps)
    plans = ligation_series(vector_length=result.backbone_length, insert_length=result.insert_length)
    text = protocol.cloning_worksheet(result, vector_name='reference', primer_set=primers,
        ligation_plans=plans, preflight=cloning_preflight(result),
        provenance={'engine_version': 'test', 'output_sha256': 'output-hash', 'parameters_sha256': 'parameter-hash'})
    assert f'{result.length} bp' in text
    assert 'output-hash' in text and 'parameter-hash' in text
    assert ('COVERAGE WARNING' in text) == bool(gaps)
    for primer in primers.primers:
        assert primer.name in text and primer.sequence in text
    for plan in plans:
        assert f'{plan.ratio:g}:1' in text


def test_protocol_exposes_duplex_errors(monkeypatch):
    plan = esd.design_extended_sequence('GCT' * 8)
    monkeypatch.setattr(protocol, 'construct_duplex', lambda _: DuplexView(top='AC', bottom='TT'))
    assert '1 positions do not pair' in protocol.bench_protocol(plan)


def test_linear_primer_coordinates_are_bounded():
    dna = clean_filler(900)
    result = design_sequencing_primers(dna, target_start=100, target_end=700, circular=False)
    assert result.primers
    assert all(0 <= p.reads_from <= len(dna) and 0 <= p.reads_to <= len(dna) for p in result.primers)
    edge = design_sequencing_primers('N' * 50 + dna, target_start=0, target_end=100, circular=False)
    assert all('N' not in p.sequence for p in edge.primers)


def test_empty_assembly_and_fragment_temperature():
    plan = esd.design_extended_sequence('GCT' * 8)
    fragment = plan.fragments[0]
    assert fragment.reverse_tm == round(melting_temperature(fragment.reverse, conditions=ANNEALING), 1)
    plan.fragments = []
    assert plan.terminal_ends == (('', 'blunt'), ('', 'blunt'))


def test_duplicate_junctions_are_rejected():
    plan = esd.design_extended_sequence(clean_filler(240), target_oligo_length=60)
    plan.fragments[1] = replace(plan.fragments[1], right_overhang=plan.fragments[0].right_overhang)
    assert any('share the overhang' in problem for problem in plan.verify())


def test_fragment_validation_failure_is_not_returned(monkeypatch):
    monkeypatch.setattr(esd.ESDResult, 'verify', lambda _: ['invalid assembly'])
    with pytest.raises(ValueError, match='invalid assembly'):
        esd.design_extended_sequence('GCT' * 8)


def test_long_single_fragment_has_length_warning():
    plan = esd.design_extended_sequence(clean_filler(210), target_oligo_length=250)
    assert plan.longest_oligo > 200
    assert any('longest oligo' in note for note in plan.warnings)
