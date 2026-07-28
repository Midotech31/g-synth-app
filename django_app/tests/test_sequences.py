"""Tests for FASTA, GenBank and SnapGene upload parsing."""
import io

import pytest
from django.core.files.uploadedfile import SimpleUploadedFile
from django.urls import reverse

from apps.sequences.parsing import ParseError, detect_format, parse_sequence_file

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
        # GenBank 1..20 is 1-based inclusive → 0-based half-open [0, 20)
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


class TestFastaParsing:
    def test_extracts_sequence(self):
        rec = parse_sequence_file(FASTA, "pFASTA.fasta")
        assert rec.name == "pFASTA"
        assert rec.length == 120
        assert rec.topology == "linear"      # FASTA carries no topology
        assert rec.annotations == []
        assert "some description" in rec.description


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
    """A vector arrives from the supplier as .dna. Asking someone to convert
    it first is asking them to use another program to use this one."""

    PET21A = "/root/.claude/uploads/7c8c4b0e-e761-5671-9d9c-44f3743c8a03/dbd85f04-pET21a.dna"

    def payload(self, name="pET-21a.dna"):
        import pathlib

        raw = pathlib.Path(self.PET21A).read_bytes()
        return SimpleUploadedFile(name, raw, content_type="application/octet-stream")

    def test_a_snapgene_file_is_recognised_by_its_bytes(self):
        """These arrive renamed as often as not."""
        import pathlib

        from apps.sequences.parsing import detect_format

        raw = pathlib.Path(self.PET21A).read_bytes()
        assert detect_format(raw, "pET-21a.dna") == "snapgene"
        assert detect_format(raw, "whatever.txt") == "snapgene"

    def test_it_parses_with_its_features(self, auth_client):
        """The reason to accept the format at all: GenBank-grade annotation."""
        response = auth_client.post(
            reverse("sequence-parse"), {"file": self.payload()}, format="multipart",
        )
        assert response.status_code == 200, response.data
        assert response.data["length"] == 5443
        assert response.data["topology"] == "circular"

        names = {a["name"] for a in response.data["annotations"]}
        assert {"AmpR", "lacI", "T7 promoter", "6xHis"} <= names

    def test_features_survive_for_every_format_that_has_them(self, auth_client):
        """The guard used to test the format name rather than the record, so
        a SnapGene file came back with its features silently dropped."""
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
