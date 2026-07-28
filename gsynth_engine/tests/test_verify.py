"""Tests for sequencing verification.

Built on a real construct: a designed insert cloned into pET-21a(+), then
reads simulated from the resulting plasmid the way a sequencing facility
would return them — in either orientation, with noise at both ends.
"""
from __future__ import annotations

import pytest

from gsynth_engine import vectors
from gsynth_engine.cloning import clone
from gsynth_engine.sequence import SequenceError, reverse_complement
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.tests.test_cloning import clean_filler
from gsynth_engine.verify import verify, verify_read


@pytest.fixture(scope="module")
def construct():
    record = vectors.sequence_of("pET-21a")
    design = design_small_sequence(clean_filler(240, 11), enzyme_pair="NdeI / XhoI")
    return clone(
        record["sequence"], design.forward, insert_reverse=design.reverse,
        left_enzyme="NdeI", right_enzyme="XhoI",
        orf_start=design.orf_start, name="pGS",
    ), design


def sanger(plasmid: str, start: int, end: int, *, noise: int = 30) -> str:
    """A read as a facility returns it: junk, sequence, junk."""
    length = len(plasmid)
    body = "".join(plasmid[i % length] for i in range(start, end))
    return "A" * noise + body + "T" * noise


class TestPlacing:
    def test_a_clean_read_matches_perfectly(self, construct):
        result, _ = construct
        read = sanger(result.plasmid, result.insert_start - 150, result.insert_end + 150)
        aligned = verify_read(result.plasmid, read, circular=True)

        assert aligned.identity == 100.0
        assert aligned.is_clean
        assert not aligned.reverse_complemented

    def test_a_reversed_read_is_recognised_and_handled(self, construct):
        """Half of all Sanger reads come back on the other strand."""
        result, _ = construct
        read = sanger(result.plasmid, result.insert_start - 150, result.insert_end + 150)
        aligned = verify_read(result.plasmid, reverse_complement(read), circular=True)

        assert aligned.reverse_complemented
        assert aligned.identity == 100.0
        assert aligned.is_clean

    def test_the_noisy_ends_are_trimmed_before_comparing(self, construct):
        """Left in, they produce differences indistinguishable from mutations."""
        result, _ = construct
        read = sanger(result.plasmid, result.insert_start, result.insert_end, noise=40)

        untrimmed = verify_read(result.plasmid, read, circular=True, trim=0)
        trimmed = verify_read(result.plasmid, read, circular=True, trim=45)

        assert untrimmed.differences
        assert trimmed.is_clean

    def test_a_read_across_the_origin_is_continuous(self, construct):
        """Split in the coordinates, continuous in the molecule."""
        result, _ = construct
        read = sanger(result.plasmid, len(result.plasmid) - 300, len(result.plasmid) + 300)
        aligned = verify_read(result.plasmid, read, circular=True)
        assert aligned.identity == 100.0

    def test_a_foreign_read_is_refused_with_a_reason(self, construct):
        result, _ = construct
        with pytest.raises(SequenceError, match="does not match the design"):
            verify_read(result.plasmid, clean_filler(400, 77), circular=True)

    def test_a_read_too_short_to_place_says_so(self, construct):
        result, _ = construct
        with pytest.raises(SequenceError, match="not enough left to place it"):
            verify_read(result.plasmid, "ACGTACGTACGT", circular=True, trim=30)

    def test_an_empty_read_is_refused(self, construct):
        result, _ = construct
        with pytest.raises(SequenceError, match="is empty"):
            verify_read(result.plasmid, "", circular=True)


class TestDifferences:
    def mutate(self, plasmid: str, position: int, base: str) -> str:
        return plasmid[:position] + base + plasmid[position + 1:]

    def test_a_substitution_is_found_at_the_right_position(self, construct):
        result, _ = construct
        at = result.insert_start + 60
        replacement = "G" if result.plasmid[at] != "G" else "C"
        mutated = self.mutate(result.plasmid, at, replacement)

        aligned = verify_read(
            result.plasmid,
            sanger(mutated, result.insert_start - 100, result.insert_end),
            circular=True,
        )
        assert len(aligned.differences) == 1
        difference = aligned.differences[0]
        assert difference.kind == "substitution"
        assert difference.position == at
        assert difference.found == replacement

    def test_a_substitution_is_reported_as_its_effect_on_the_protein(self, construct):
        """What the user needs to know is whether it matters."""
        result, design = construct
        at = result.insert_start + design.orf_start + 30
        replacement = "G" if result.plasmid[at] != "G" else "C"
        mutated = self.mutate(result.plasmid, at, replacement)

        aligned = verify_read(
            result.plasmid,
            sanger(mutated, result.insert_start - 100, result.insert_end),
            circular=True,
            coding_start=result.insert_start + design.orf_start,
            coding_end=result.insert_end,
        )
        difference = aligned.differences[0]
        assert difference.residue == 11          # base 30 of the frame
        assert difference.silent in (True, False)
        assert difference.from_residue
        assert "residue" in difference.description

    def test_a_deletion_is_found(self, construct):
        result, _ = construct
        at = result.insert_start + 80
        deleted = result.plasmid[:at] + result.plasmid[at + 1:]

        aligned = verify_read(
            result.plasmid,
            sanger(deleted, result.insert_start - 100, result.insert_end - 1),
            circular=True,
        )
        assert any(d.kind == "deletion" for d in aligned.differences)

    def test_an_insertion_is_found(self, construct):
        result, _ = construct
        at = result.insert_start + 80
        inserted = result.plasmid[:at] + "G" + result.plasmid[at:]

        aligned = verify_read(
            result.plasmid,
            sanger(inserted, result.insert_start - 100, result.insert_end),
            circular=True,
        )
        assert any(d.kind == "insertion" for d in aligned.differences)

    def test_a_poor_read_is_flagged_rather_than_believed(self, construct):
        """A dozen differences in one read is a bad trace, not a dozen mutations."""
        result, _ = construct
        body = list(result.plasmid[result.insert_start : result.insert_end])
        for i in range(20, len(body), 30):
            body[i] = "G" if body[i] != "G" else "C"

        aligned = verify_read(
            result.plasmid, "A" * 30 + "".join(body) + "T" * 30,
            circular=True,
        )
        assert len(aligned.differences) > 3
        assert aligned.identity < 98
        assert any("trace quality" in note for note in aligned.warnings)

    def test_a_hopeless_read_is_refused_rather_than_misreported(self, construct):
        """At one error in nine there is nothing to anchor on. Saying so beats
        reporting forty mutations that are really a failed run."""
        result, _ = construct
        body = list(result.plasmid[result.insert_start : result.insert_end])
        for i in range(5, len(body), 9):
            body[i] = "G" if body[i] != "G" else "C"

        with pytest.raises(SequenceError, match="does not match the design"):
            verify_read(
                result.plasmid, "A" * 30 + "".join(body) + "T" * 30, circular=True,
            )


class TestReport:
    def test_two_reads_covering_the_insert_verify_it(self, construct):
        result, _ = construct
        reads = {
            "T7-F": sanger(result.plasmid, result.insert_start - 200, result.insert_start + 300),
            "T7-R": reverse_complement(
                sanger(result.plasmid, result.insert_start + 200, result.insert_end + 200)
            ),
        }
        report = verify(
            result.plasmid, reads, circular=True,
            region=(result.insert_start, result.insert_end),
        )
        assert report.is_verified
        assert report.fully_covered
        assert report.coverage == 100.0
        assert len(report.reads) == 2

    def test_a_gap_in_the_coverage_is_reported(self, construct):
        """Half the insert read is not the insert verified."""
        result, _ = construct
        report = verify(
            result.plasmid,
            {"one": sanger(result.plasmid, result.insert_start - 100, result.insert_start + 60)},
            circular=True, region=(result.insert_start, result.insert_end),
        )
        assert not report.fully_covered
        assert report.coverage < 100

    def test_one_bad_read_does_not_stop_the_others(self, construct):
        result, _ = construct
        report = verify(
            result.plasmid,
            {
                "good": sanger(result.plasmid, result.insert_start - 100, result.insert_end),
                "junk": clean_filler(400, 91),
            },
            circular=True, region=(result.insert_start, result.insert_end),
        )
        assert len(report.reads) == 1
        assert report.warnings
        assert report.is_verified

    def test_no_reads_means_not_verified(self, construct):
        result, _ = construct
        report = verify(result.plasmid, [], circular=True)
        assert not report.is_verified

    def test_the_same_change_seen_twice_is_listed_once(self, construct):
        result, _ = construct
        at = result.insert_start + 60
        mutated = result.plasmid[:at] + ("G" if result.plasmid[at] != "G" else "C") + result.plasmid[at + 1:]
        window = (result.insert_start - 100, result.insert_end)

        report = verify(
            result.plasmid,
            {"a": sanger(mutated, *window), "b": sanger(mutated, *window)},
            circular=True, region=(result.insert_start, result.insert_end),
        )
        assert len(report.differences) == 1
        assert not report.is_verified
