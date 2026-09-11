"""Independent sequence extraction checks for coding feature transport."""
import io

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from gsynth_engine.cloning import _flip_annotations, _remap_annotations, linearise
from gsynth_engine.genbank import to_genbank


@pytest.mark.parametrize('direction', [-1, 1])
@pytest.mark.parametrize('wrapped', [False, True])
def test_export_uses_the_actual_coding_interval(direction, wrapped):
    dna = 'ACGT' * 15
    start, end = (52, 70) if wrapped else (2, 20)
    feature = {'name': 'coding interval', 'type': 'CDS', 'start': start, 'end': end,
               'translation_start': start + 1, 'translation_end': end - 2, 'direction': direction}
    record = SeqIO.read(io.StringIO(to_genbank(dna, features=[feature], circular=True)), 'genbank')
    cds = next(f for f in record.features if f.type == 'CDS')
    expected = Seq((dna * 2)[start + 1:end - 2])
    if direction == -1:
        expected = expected.reverse_complement()
    assert str(cds.extract(record.seq)).upper() == str(expected)


@pytest.mark.parametrize('direction', [-1, 0, 1])
def test_mirroring_twice_preserves_wrapped_coding_bounds_and_strand(direction):
    feature = {'start': 52, 'end': 70, 'translation_start': 53, 'translation_end': 68, 'direction': direction}
    flipped = _flip_annotations([feature], 60)
    assert 0 <= flipped[0]['start'] < 60
    assert _flip_annotations(flipped, 60) == [feature]


def test_backbone_annotation_transport_moves_coding_bounds_with_the_feature():
    vector = 'A' * 40 + 'CATATG' + 'A' * 15 + 'CTCGAG' + 'ACGTTGCAAGCTTAGCGATCCGTACAGTTCGAGCGTACCATGGCAGTACGATCAGTGCAT'
    backbone = linearise(vector, left_enzyme='NdeI', right_enzyme='XhoI')
    feature = {'name': 'coding', 'type': 'CDS', 'start': 80, 'end': 98,
               'translation_start': 81, 'translation_end': 96, 'direction': -1}
    moved = _remap_annotations([feature], backbone, backbone.length + 30)[0]
    expected = Seq(vector[81:96]).reverse_complement()
    actual = Seq(backbone.top[moved['translation_start']:moved['translation_end']]).reverse_complement()
    assert actual == expected
