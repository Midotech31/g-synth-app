import io

import pytest
import sbol3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from apps.sequences.parsing import ParseError, parse_sequence_file
from apps.sequences.sbol import read_sbol3, to_sbol3


@pytest.mark.parametrize('direction', [-1, 0, 1])
@pytest.mark.parametrize('wrapped', [False, True])
def test_sbol_preserves_annotation_evidence_and_coding_bounds(direction, wrapped):
    start, end = (52, 70) if wrapped else (2, 20)
    annotation = {'name': 'Candidate interval', 'type': 'CDS', 'start': start, 'end': end,
        'direction': direction, 'color': '#123456', 'inferred': True,
        'basis': 'Candidate from source record; function unconfirmed.',
        'regulatory_class': '', 'translation_start': start + 1, 'translation_end': end - 2}
    text = to_sbol3('ACGT' * 15, features=[annotation], circular=wrapped)
    record = parse_sequence_file(text, 'test.sbol.json')
    actual = record.annotations[0]
    for key, value in annotation.items():
        assert getattr(actual, key) == value


def test_sbol_rejects_discontinuous_feature_instead_of_flattening():
    text = to_sbol3('ACGT' * 15, features=[{'name': 'split', 'start': 2, 'end': 10}])
    document = read_sbol3(text, 'sbol3-jsonld')
    component = next(item for item in document if isinstance(item, sbol3.Component))
    component.features[0].locations.append(sbol3.Range(component.sequences[0], 30, 40))
    with pytest.raises(ParseError, match='Discontinuous'):
        parse_sequence_file(document.write_string(sbol3.JSONLD), 'test.sbol.json')


@pytest.mark.parametrize('feature,circular', [({'start': -1, 'end': 3}, False),
    ({'start': 2, 'end': 10}, False), ({'start': 2, 'end': 11}, True)])
def test_sbol_export_rejects_invalid_coordinates(feature, circular):
    with pytest.raises(ValueError):
        to_sbol3('ACGTACGT', features=[feature], circular=circular)


def test_sbol_rejects_invalid_coding_bounds():
    text = to_sbol3('ACGT' * 15, features=[{'start': 2, 'end': 20, 'type': 'CDS',
        'translation_start': 1, 'translation_end': 10}])
    with pytest.raises(ParseError, match='coding bounds'):
        parse_sequence_file(text, 'test.sbol.json')


@pytest.mark.parametrize('qualifiers', [{'transl_table': ['2']}, {'transl_except': ['(pos:4..6,aa:Sec)']}])
def test_genbank_unsupported_translation_is_explicit(qualifiers):
    record = SeqRecord(Seq('ATGTGATAA'), id='test', annotations={'molecule_type': 'DNA'})
    record.features = [SeqFeature(SimpleLocation(0, 9, strand=1), type='CDS', qualifiers=qualifiers)]
    out = io.StringIO()
    SeqIO.write(record, out, 'genbank')
    with pytest.raises(ParseError, match='unsupported genetic code'):
        parse_sequence_file(out.getvalue(), 'test.gb')
