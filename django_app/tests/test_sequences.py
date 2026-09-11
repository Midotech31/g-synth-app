import io
import json
import struct
from pathlib import Path
from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from django.core.files.uploadedfile import SimpleUploadedFile
from django.urls import reverse

from apps.sequences.parsing import (
    ParseError,
    _annotations_from,
    detect_format,
    parse_sequence_file,
)
from apps.sequences.sbol import to_sbol3


@pytest.mark.parametrize('strand', [1, -1])
@pytest.mark.parametrize('offset', [0, 1, 2])
@pytest.mark.parametrize('wrapped', [False, True])
def test_cds_frame_survives_import_export_against_biopython(strand, offset, wrapped):
    from Bio import SeqIO
    from Bio.SeqFeature import CompoundLocation

    from gsynth_engine.genbank import to_genbank
    sequence = Seq('ACGT' * 15)
    record = SeqRecord(sequence, id='frame', annotations={'molecule_type': 'DNA', 'topology': 'circular'})
    parts = [SimpleLocation(50, 60, strand=strand), SimpleLocation(0, 10, strand=strand)]
    location = CompoundLocation(parts if strand == 1 else parts[::-1]) if wrapped else SimpleLocation(5, 25, strand=strand)
    original = SeqFeature(location, type='CDS', qualifiers={'codon_start': [str(offset + 1)]})
    record.features = [original]
    expected = str(original.extract(sequence))[offset:]
    text = io.StringIO()
    SeqIO.write(record, text, 'genbank')
    imported = parse_sequence_file(text.getvalue().encode(), 'frame.gb')
    exported = to_genbank(imported.sequence, circular=True, features=imported.to_dict()['annotations'])
    independent = SeqIO.read(io.StringIO(exported), 'genbank')
    cds = next(f for f in independent.features if f.type == 'CDS')
    assert str(cds.extract(independent.seq)).upper() == expected
    assert cds.qualifiers['codon_start'] == ['1']


def test_discontinuous_feature_is_not_silently_flattened():
    from Bio.SeqFeature import CompoundLocation
    record = SeqRecord(Seq('ACGT' * 15), id='split')
    record.features = [SeqFeature(CompoundLocation([SimpleLocation(1, 5), SimpleLocation(10, 15)]), type='CDS')]
    with pytest.raises(ParseError, match='discontinuous'):
        _annotations_from(record)

GENBANK = """LOCUS       pDEMO                    120 bp    DNA     circular SYN 01-JAN-2026
DEFINITION  Demonstration plasmid for the G-Synth viewer.
ACCESSION   pDEMO
FEATURES             Location/Qualifiers
     source          1..120
                     /organism="synthetic construct"
     promoter        1..20
                     /label="T7 promoter"
     CDS             21..80
                     /label="GFP fragment"
                     /codon_start=1
     terminator      complement(90..115)
                     /label="rrnB T1"
ORIGIN
        1 taatacgact cactataggg atggtgagca agggcgagga gctgttcacc ggggtggtgc
       61 ccatcctggt cgagctggac ggcgacgtaa acggccacaa gttcagcgtg tccggcgagg
//
"""

FASTA = """>pFASTA some description here
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC
GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC
"""


def _snapgene_pet21a() -> bytes:

    path = Path(__file__).resolve().parents[2] / "gsynth_engine/vector_data/pET-21a.json"
    record = json.loads(path.read_text(encoding="utf-8"))

    def packet(kind: int, data: bytes) -> bytes:
        return bytes([kind]) + struct.pack(">I", len(data)) + data

    cookie = packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 1, 1))
    sequence = packet(0x00, b"\x01" + record["sequence"].encode("ascii"))

    features = ElementTree.Element("Features")
    for annotation in record["annotations"]:
        feature = ElementTree.SubElement(
            features,
            "Feature",
            name=annotation["name"],
            type=annotation["type"],
            directionality="2" if annotation["direction"] == -1 else "1",
        )
        ElementTree.SubElement(
            feature,
            "Segment",
            range=f"{annotation['start'] + 1}-{annotation['end']}",
            type="standard",
        )
    return cookie + sequence + packet(0x0A, ElementTree.tostring(features))


class TestFormatDetection:
    def test_detects_genbank_by_content(self):
        assert detect_format(GENBANK) == "genbank"

    def test_detects_fasta_by_content(self):
        assert detect_format(FASTA) == "fasta"

    def test_falls_back_to_extension(self):
        assert detect_format("ACGTACGT\n", "mystery.fasta") == "fasta"

    def test_rejects_unknown(self):
        with pytest.raises(ParseError):
            detect_format("just some prose", "notes.txt")


class TestGenBankParsing:
    def test_extracts_sequence_and_metadata(self):
        rec = parse_sequence_file(GENBANK, "pDEMO.gb")
        assert rec.name == "pDEMO"
        assert rec.length == 120
        assert len(rec.sequence) == 120
        assert rec.topology == "circular"
        assert rec.source_format == "genbank"
        assert "Demonstration plasmid" in rec.description

    def test_gc_content_matches_manual_count(self):
        rec = parse_sequence_file(GENBANK, "pDEMO.gb")
        gc = sum(1 for b in rec.sequence if b in "GC")
        assert rec.gc_content == pytest.approx(100.0 * gc / len(rec.sequence), abs=0.05)

    def test_annotations_are_viewer_ready(self):
        rec = parse_sequence_file(GENBANK, "pDEMO.gb")
        by_name = {a.name: a for a in rec.annotations}
        assert "T7 promoter" in by_name
        assert "GFP fragment" in by_name
        assert "rrnB T1" in by_name

        promoter = by_name["T7 promoter"]

        assert (promoter.start, promoter.end) == (0, 20)
        assert promoter.direction == 1
        assert promoter.color.startswith("#")

        terminator = by_name["rrnB T1"]
        assert terminator.direction == -1, "complement(...) must read as reverse strand"

    def test_source_feature_is_not_drawn(self):
        rec = parse_sequence_file(GENBANK, "pDEMO.gb")
        assert all(a.type != "source" for a in rec.annotations)

    def test_annotations_sorted_by_start(self):
        rec = parse_sequence_file(GENBANK, "pDEMO.gb")
        starts = [a.start for a in rec.annotations]
        assert starts == sorted(starts)

    def test_unstranded_feature_stays_unstranded(self):
        record = SeqRecord(Seq("ACGT"))
        record.features = [
            SeqFeature(SimpleLocation(0, 4, strand=None), type="misc_feature")
        ]

        assert _annotations_from(record)[0].direction == 0

    def test_a_feature_crossing_the_circular_origin_keeps_its_true_span(self):
        text = GENBANK.replace(
            "     promoter        1..20\n",
            "     promoter        join(111..120,1..10)\n",
        )
        record = parse_sequence_file(text, "wrapped.gb")
        feature = next(a for a in record.annotations if a.name == "T7 promoter")
        assert (feature.start, feature.end) == (110, 130)
        assert feature.end - feature.start == 20


class TestFastaParsing:
    def test_extracts_sequence(self):
        rec = parse_sequence_file(FASTA, "pFASTA.fasta")
        assert rec.name == "pFASTA"
        assert rec.length == 120
        assert rec.topology == "linear"
        assert rec.annotations == []
        assert "some description" in rec.description


class TestSbol3Parsing:
    def test_validated_jsonld_round_trip_preserves_sequence_topology_and_features(self):
        document = to_sbol3(
            "ATGGTGAGCAAGGGCGAGGA",
            name="pGS test",
            description="SBOL interchange construct",
            circular=True,
            features=[{
                "name": "coding region",
                "type": "CDS",
                "start": 0,
                "end": 18,
                "direction": 1,
            }],
        )

        record = parse_sequence_file(document, "pGS_test.sbol.json")

        assert record.source_format == "sbol3"
        assert record.name == "pGS test"
        assert record.topology == "circular"
        assert record.annotations[0].type == "CDS"
        assert (record.annotations[0].start, record.annotations[0].end) == (0, 18)

    def test_detects_turtle_and_jsonld_by_content(self):
        document = to_sbol3("ATGC", name="test")
        assert detect_format(document, "renamed.txt") == "sbol3-jsonld"
        assert detect_format(
            "@prefix sbol: <http://sbols.org/v3#> .", "renamed.txt"
        ) == "sbol3-turtle"


class TestParseFailures:
    def test_empty_file(self):
        with pytest.raises(ParseError, match="empty"):
            parse_sequence_file("   ")

    def test_garbage_file(self):
        with pytest.raises(ParseError):
            parse_sequence_file("this is not a sequence file", "notes.txt")

    def test_oversized_file(self):
        with pytest.raises(ParseError, match="larger than"):
            parse_sequence_file(b"A" * (11 * 1024 * 1024), "huge.fasta")


@pytest.mark.django_db
class TestParseEndpoint:
    url = "/api/sequences/parse/"

    def test_requires_auth(self, api_client):
        assert api_client.post(self.url).status_code == 401

    def test_uploads_genbank(self, auth_client):
        upload = io.BytesIO(GENBANK.encode())
        upload.name = "pDEMO.gb"
        r = auth_client.post(self.url, {"file": upload}, format="multipart")
        assert r.status_code == 200, r.data
        assert r.data["name"] == "pDEMO"
        assert r.data["topology"] == "circular"
        assert len(r.data["annotations"]) == 3
        assert {"name", "start", "end", "direction", "color"} <= set(r.data["annotations"][0])

    def test_uploads_fasta(self, auth_client):
        upload = io.BytesIO(FASTA.encode())
        upload.name = "pFASTA.fasta"
        r = auth_client.post(self.url, {"file": upload}, format="multipart")
        assert r.status_code == 200
        assert r.data["length"] == 120

    def test_uploads_sbol3_with_annotations(self, auth_client):
        text = to_sbol3(
            "ATGGTGAGCAAG",
            name="pSBOL",
            features=[{"name": "gene", "type": "gene", "start": 0, "end": 12}],
        )
        upload = io.BytesIO(text.encode())
        upload.name = "pSBOL.sbol.json"
        response = auth_client.post(self.url, {"file": upload}, format="multipart")

        assert response.status_code == 200, response.data
        assert response.data["source_format"] == "sbol3"
        assert response.data["annotations"][0]["name"] == "gene"

    def test_missing_file_is_a_clear_error(self, auth_client):
        r = auth_client.post(self.url, {}, format="multipart")
        assert r.status_code == 400
        assert "file" in r.data["detail"].lower()

    def test_garbage_upload_reports_why(self, auth_client):
        upload = io.BytesIO(b"definitely not a sequence")
        upload.name = "notes.txt"
        r = auth_client.post(self.url, {"file": upload}, format="multipart")
        assert r.status_code == 400
        assert "FASTA" in r.data["detail"]


@pytest.mark.django_db
class TestSnapGeneImport:


    def payload(self, name="pET-21a.dna"):
        return SimpleUploadedFile(
            name, _snapgene_pet21a(), content_type="application/octet-stream"
        )

    def test_a_snapgene_file_is_recognised_by_its_bytes(self):

        raw = _snapgene_pet21a()
        assert detect_format(raw, "pET-21a.dna") == "snapgene"
        assert detect_format(raw, "whatever.txt") == "snapgene"

    def test_it_parses_with_its_features(self, auth_client):

        response = auth_client.post(
            reverse("sequence-parse"), {"file": self.payload()}, format="multipart",
        )
        assert response.status_code == 200, response.data
        assert response.data["length"] == 5443
        assert response.data["topology"] == "circular"

        names = {a["name"] for a in response.data["annotations"]}
        assert {"AmpR", "lacI", "T7 promoter", "6xHis"} <= names

    def test_features_survive_for_every_format_that_has_them(self, auth_client):

        response = auth_client.post(
            reverse("sequence-parse"), {"file": self.payload()}, format="multipart",
        )
        assert len(response.data["annotations"]) == 14

    def test_fasta_still_has_none_and_still_parses(self, auth_client):
        upload = SimpleUploadedFile("x.fasta", b">x desc\nACGTACGTAA\n", content_type="text/plain")
        response = auth_client.post(
            reverse("sequence-parse"), {"file": upload}, format="multipart",
        )
        assert response.status_code == 200
        assert response.data["annotations"] == []

    def test_a_binary_file_that_is_not_snapgene_is_refused_clearly(self, auth_client):
        upload = SimpleUploadedFile("junk.bin", bytes(range(200)), content_type="application/octet-stream")
        response = auth_client.post(
            reverse("sequence-parse"), {"file": upload}, format="multipart",
        )
        assert response.status_code == 400
        assert "not a SnapGene file" in response.data["detail"]
