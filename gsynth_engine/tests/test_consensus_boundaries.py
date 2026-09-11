import importlib
from dataclasses import replace

import pytest

from gsynth_engine.chromatogram import Chromatogram
from gsynth_engine.sequence import SequenceError
from gsynth_engine.tests.test_cloning import clean_filler

verification = importlib.import_module('gsynth_engine.verify')


def trace(sequence):
    return Chromatogram(sequence, [40] * len(sequence), [], {})


@pytest.mark.parametrize('design,region,message', [('', None, 'empty'), ('ACGT', (0, 5), 'outside')])
def test_consensus_rejects_invalid_reference(design, region, message):
    with pytest.raises(SequenceError, match=message):
        verification.assemble_consensus(design, {}, region=region)


def test_empty_and_unplaced_reads_do_not_create_coverage():
    dna = clean_filler(200)
    result = verification.assemble_consensus(dna, {
        'empty': trace(''), 'short': trace('ACGT'), 'unplaced': trace('A' * 100),
    })
    assert result.coverage == 0
    assert result.sequence == 'N' * len(dna)
    assert result.gaps == [(0, len(dna))]
    assert len(result.warnings) == 3


def test_trimmed_partial_consensus_preserves_missing_intervals():
    dna = clean_filler(200)
    read = replace(trace(dna[40:160]), quality=[0] * 10 + [40] * 100 + [0] * 10)
    result = verification.assemble_consensus(dna, {'read': read}, trim_quality=13)
    assert result.sequence[50:150] == dna[50:150]
    assert result.gaps == [(0, 50), (150, 200)]
    assert result.coverage == 50


@pytest.mark.parametrize('change', ['insertion', 'deletion'])
def test_consensus_handles_indels_and_region_limits(change):
    dna = clean_filler(240)
    read = dna[:120] + ('A' + dna[120:] if change == 'insertion' else dna[121:])
    result = verification.assemble_consensus(dna, {'read': trace(read)}, region=(60, 180))
    assert len(result.sequence) == 120
    assert result.sequence[:50] == dna[60:110]
    assert all(60 <= position.position < 180 for position in result.positions)
    if change == 'deletion':
        assert result.gaps


def test_minimal_overlap_is_not_admitted(monkeypatch):
    dna = clean_filler(100)
    monkeypatch.setattr(verification, '_locate', lambda *args, **kwargs: (95, False))
    with pytest.raises(SequenceError, match='overlaps only 5'):
        verification.verify_read(dna, dna, trim=0)
    result = verification.assemble_consensus(dna, {'read': trace(dna)})
    assert result.coverage == 0
    assert any('overlaps only 5' in warning for warning in result.warnings)


def test_coding_effect_outside_or_incomplete_codon():
    assert verification._coding_effect('ATGA', 0, 'C', 1, 4) == {}
    assert verification._coding_effect('ATGA', 3, 'C', 0, 4) == {}
    with pytest.raises(SequenceError, match='empty'):
        verification.verify_read('', 'ATG', trim=0)


def test_incomplete_alignment_is_rejected():
    with pytest.raises(SequenceError, match='cannot be verified completely'):
        verification._align('CGCAACGCGTGAATTAT', 'ACTTTATGCGTTATAAA', 5, circular=False)


def test_consensus_keeps_valid_reads_when_one_alignment_fails(monkeypatch):
    dna = clean_filler(200)
    original = verification._align
    def align(design, read, *args, **kwargs):
        if read != design:
            raise SequenceError('Incomplete alignment')
        return original(design, read, *args, **kwargs)
    monkeypatch.setattr(verification, '_align', align)
    result = verification.assemble_consensus(dna, {'valid': trace(dna), 'invalid': trace(dna[20:180])})
    assert result.sequence == dna
    assert result.warnings == ['invalid: Incomplete alignment']
