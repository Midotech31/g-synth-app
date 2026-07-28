"""Tests for GenBank and FASTA output.

Writing a format is easy to do almost right, and almost right is worse than
wrong: a LOCUS line whose columns are off by two produces a record that
parses without complaint and comes back linear. So the tests here do not
check that the text looks plausible — they parse it back with Biopython, an
implementation that had no part in writing it, and compare.
"""
from __future__ import annotations

import io

import pytest
from Bio import SeqIO

from gsynth_engine import vectors
from gsynth_engine.cloning import clone
from gsynth_engine.genbank import (
    Feature,
    oligos_to_fasta,
    to_fasta,
    to_genbank,
)
from gsynth_engine.ssd import design_small_sequence

SEQUENCE = "ATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCTTAA"
FEATURES = [
    {"name": "6xHis", "type": "CDS", "start": 12, "end": 30, "direction": 1},
    {"name": "thrombin", "type": "misc_feature", "start": 39, "end": 57, "direction": -1},
]


def parse(text: str):
    """Read a record back with an independent parser."""
    return SeqIO.read(io.StringIO(text), "genbank")


# ── The round trip ──────────────────────────────────────────────────────────


class TestRoundTrip:
    def test_the_sequence_survives(self):
        record = parse(to_genbank(SEQUENCE, name="cassette"))
        assert str(record.seq).upper() == SEQUENCE

    def test_the_name_survives(self):
        record = parse(to_genbank(SEQUENCE, name="pGS-EntA"))
        assert record.name == "pGS-EntA"

    def test_circular_topology_survives(self):
        """The failure this file exists to prevent: a plasmid read as linear."""
        record = parse(to_genbank(SEQUENCE, name="plasmid", circular=True))
        assert record.annotations["topology"] == "circular"

    def test_linear_stays_linear(self):
        record = parse(to_genbank(SEQUENCE, name="fragment", circular=False))
        assert record.annotations["topology"] == "linear"

    def test_features_survive_with_their_coordinates(self):
        record = parse(to_genbank(SEQUENCE, features=FEATURES))
        drawn = [f for f in record.features if f.type != "source"]
        assert len(drawn) == len(FEATURES)

        for written, read in zip(FEATURES, drawn):
            assert int(read.location.start) == written["start"]
            assert int(read.location.end) == written["end"]

    def test_features_keep_their_strand(self):
        """A resistance gene on the wrong strand is nonsense on the map."""
        record = parse(to_genbank(SEQUENCE, features=FEATURES))
        drawn = [f for f in record.features if f.type != "source"]
        assert drawn[0].location.strand == 1
        assert drawn[1].location.strand == -1

    def test_feature_labels_survive(self):
        record = parse(to_genbank(SEQUENCE, features=FEATURES))
        labels = [
            f.qualifiers.get("label", [""])[0]
            for f in record.features if f.type != "source"
        ]
        assert labels == ["6xHis", "thrombin"]

    def test_feature_types_survive(self):
        record = parse(to_genbank(SEQUENCE, features=FEATURES))
        types = [f.type for f in record.features if f.type != "source"]
        assert types == ["CDS", "misc_feature"]

    def test_the_description_survives(self):
        record = parse(to_genbank(SEQUENCE, description="EntA in pET-21a(+)"))
        assert "EntA in pET-21a(+)" in record.description


# ── Real payloads ───────────────────────────────────────────────────────────


class TestRealConstructs:
    def test_a_recombinant_plasmid_round_trips(self):
        """The output that matters: what the user opens in SnapGene."""
        record = vectors.sequence_of("pET-21a")
        design = design_small_sequence(
            "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGC", enzyme_pair="NdeI / XhoI"
        )
        result = clone(
            record["sequence"], design.forward, insert_reverse=design.reverse,
            left_enzyme="NdeI", right_enzyme="XhoI",
            vector_annotations=record["annotations"], name="pGS-EntA",
        )

        text = to_genbank(
            result.plasmid, name="pGS-EntA", circular=True,
            features=result.annotations,
        )
        parsed = parse(text)

        assert str(parsed.seq).upper() == result.plasmid
        assert parsed.annotations["topology"] == "circular"

        labels = {
            f.qualifiers.get("label", [""])[0]
            for f in parsed.features if f.type != "source"
        }
        assert {"AmpR", "ori", "lacI"} <= labels

    def test_a_bundled_vector_round_trips(self):
        record = vectors.sequence_of("pET-21a")
        parsed = parse(
            to_genbank(
                record["sequence"], name=record["name"], circular=True,
                features=record["annotations"],
            )
        )
        assert len(parsed.seq) == 5443
        assert len([f for f in parsed.features if f.type != "source"]) == 14

    def test_a_five_kilobase_record_wraps_correctly(self):
        record = vectors.sequence_of("pET-21a")
        text = to_genbank(record["sequence"], name="pET21a")
        assert str(parse(text).seq).upper() == record["sequence"]


# ── Edge cases ──────────────────────────────────────────────────────────────


class TestEdges:
    def test_a_feature_across_the_origin_becomes_a_join(self):
        """On a circle a feature can end past the last base."""
        text = to_genbank(
            SEQUENCE, circular=True,
            features=[{"name": "wrap", "start": 55, "end": 65, "direction": 1}],
        )
        assert "join(" in text
        record = parse(text)
        wrapped = [f for f in record.features if f.type != "source"][0]
        assert len(wrapped.location.parts) == 2

    def test_an_empty_feature_is_skipped(self):
        """Zero-length features come from clipping and draw as artefacts."""
        record = parse(
            to_genbank(SEQUENCE, features=[{"name": "gone", "start": 10, "end": 10}])
        )
        assert [f for f in record.features if f.type != "source"] == []

    def test_a_long_name_does_not_break_the_locus_line(self):
        record = parse(
            to_genbank(SEQUENCE, name="a-very-long-construct-name-indeed")
        )
        assert str(record.seq).upper() == SEQUENCE

    def test_a_name_with_spaces_does_not_break_the_locus_line(self):
        record = parse(to_genbank(SEQUENCE, name="my construct 2"))
        assert str(record.seq).upper() == SEQUENCE

    def test_output_is_the_same_every_time(self):
        """Someone will diff two exports; a clock in the file defeats that."""
        first = to_genbank(SEQUENCE, name="x", features=FEATURES)
        second = to_genbank(SEQUENCE, name="x", features=FEATURES)
        assert first == second

    def test_a_dataclass_feature_works_as_well_as_a_dict(self):
        record = parse(
            to_genbank(
                SEQUENCE,
                features=[Feature(name="tag", type="CDS", start=0, end=9)],
            )
        )
        assert [f.type for f in record.features if f.type != "source"] == ["CDS"]

    def test_an_empty_sequence_still_produces_a_record(self):
        assert to_genbank("", name="empty").endswith("//\n")


# ── FASTA ───────────────────────────────────────────────────────────────────


class TestFasta:
    def test_it_round_trips(self):
        record = SeqIO.read(io.StringIO(to_fasta(SEQUENCE, name="cassette")), "fasta")
        assert str(record.seq) == SEQUENCE
        assert record.id == "cassette"

    def test_the_description_is_kept(self):
        text = to_fasta(SEQUENCE, name="cassette", description="EntA construct")
        record = SeqIO.read(io.StringIO(text), "fasta")
        assert record.description.endswith("EntA construct")

    def test_lines_are_wrapped(self):
        lines = to_fasta("A" * 200, width=70).splitlines()
        assert all(len(line) <= 70 for line in lines[1:])

    def test_every_oligo_becomes_one_entry(self):
        """Suppliers take a FASTA upload; retyping thirty names is where
        transcription errors come from."""
        oligos = [
            {"Name": "F1_F", "Sequence (5'->3')": "ATGCATGC"},
            {"Name": "F1_R", "Sequence (5'->3')": "GCATGCAT"},
        ]
        records = list(SeqIO.parse(io.StringIO(oligos_to_fasta(oligos)), "fasta"))
        assert [r.id for r in records] == ["F1_F", "F1_R"]
        assert str(records[0].seq) == "ATGCATGC"

    def test_an_empty_order_sheet_is_empty(self):
        assert oligos_to_fasta([]) == ""
