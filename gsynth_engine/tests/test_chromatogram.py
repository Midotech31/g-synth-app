"""Tests for reading Sanger trace files.

The point of reading the trace rather than the letters is that a difference
from the design means two completely different things depending on the peak
underneath it. A mismatch at Q50 is a mutation and someone must decide what
to do about it; the same mismatch at Q8 is the basecaller guessing between a
peak and its shoulder, and chasing it wastes a week.

These build ABIF files rather than shipping one, because a fixture from a
particular instrument tests that instrument. The writer here follows the
published format, so the reader is checked against the specification and not
against whatever one facility happened to send.
"""
import struct

import pytest

from gsynth_engine.chromatogram import (
    DEFAULT_TRIM_QUALITY,
    read_ab1,
    summarise,
)
from gsynth_engine.sequence import SequenceError

READ = "ACGTTGCAAGGCTTAGCCATGGATCCGTTAACCGGTTAAGCTTCAGGATCCAAGCTTGGCC"


def build_ab1(sequence: str, quality: list[int] | None = None, *,
              order: str = "GATC", peaks: list[int] | None = None,
              spacing: int = 12, traces: dict[str, list[int]] | None = None,
              omit: set[str] = frozenset()) -> bytes:
    """A valid ABIF file, written to the published layout.

    Header, then a directory of tagged entries, then their data. Anything of
    four bytes or fewer lives inside the entry's offset field rather than
    being pointed to by it — the case a reader is most likely to get wrong.
    """
    quality = quality if quality is not None else [40] * len(sequence)
    peaks = peaks if peaks is not None else [spacing * (i + 1) for i in range(len(sequence))]

    if traces is None:
        samples = (len(sequence) + 2) * spacing
        traces = {b: [0] * samples for b in "ACGT"}
        for i, base in enumerate(sequence):
            at = peaks[i]
            if base in traces and at < samples:
                for d in range(-3, 4):                 # a peak with shoulders
                    if 0 <= at + d < samples:
                        traces[base][at + d] += max(0, 800 - 220 * abs(d))

    entries: list[tuple[str, int, int, int, bytes]] = []

    def add(tag: str, number: int, kind: int, count: int, payload: bytes):
        if tag not in omit:
            entries.append((tag, number, kind, count, payload))

    add("PBAS", 1, 2, len(sequence), sequence.encode("ascii"))
    add("PCON", 1, 1, len(quality), bytes(min(255, max(0, q)) for q in quality))
    add("PLOC", 1, 4, len(peaks), struct.pack(f">{len(peaks)}h", *peaks))
    add("FWO_", 1, 2, 4, order.encode("ascii"))
    for i, base in enumerate(order[:4]):
        channel = traces.get(base, [])
        add("DATA", 9 + i, 4, len(channel), struct.pack(f">{len(channel)}h", *channel))

    directory_at = 128
    data_at = directory_at + len(entries) * 28

    blob = bytearray(b"\x00" * data_at)
    blob[0:4] = b"ABIF"
    blob[4:6] = struct.pack(">h", 101)

    # The header's own entry describes the directory that follows it.
    blob[6:34] = struct.pack(">4sihhii4si", b"tdir", 1, 1023, 28, len(entries),
                             len(entries) * 28, struct.pack(">i", directory_at), 0)

    for i, (tag, number, kind, count, payload) in enumerate(entries):
        if len(payload) <= 4:
            offset = payload.ljust(4, b"\x00")
        else:
            offset = struct.pack(">i", len(blob))
            blob += payload
        blob[directory_at + i * 28 : directory_at + (i + 1) * 28] = struct.pack(
            ">4sihhii4si", tag.encode("ascii"), number, kind,
            {1: 1, 2: 1, 4: 2}[kind], count, len(payload), offset, 0)

    return bytes(blob)


class TestReadingATrace:
    def test_base_calls_survive_the_round_trip(self):
        trace = read_ab1(build_ab1(READ))
        assert trace.sequence == READ

    def test_quality_is_kept_per_base(self):
        quality = [10 + (i % 40) for i in range(len(READ))]
        trace = read_ab1(build_ab1(READ, quality))
        assert trace.quality == quality
        assert trace.quality_at(3) == quality[3]

    def test_quality_past_either_end_is_zero_not_an_error(self):
        """Callers index this from alignment coordinates, which can run off."""
        trace = read_ab1(build_ab1(READ))
        assert trace.quality_at(-1) == 0
        assert trace.quality_at(10_000) == 0

    def test_all_four_channels_are_read(self):
        trace = read_ab1(build_ab1(READ))
        assert set(trace.traces) == set("ACGT")
        assert all(trace.traces[b] for b in "ACGT")
        assert trace.sample_count == len(trace.traces["A"])

    def test_channels_follow_the_files_own_base_order(self):
        """FWO_ says which base DATA9-12 carry, and it is not always GATC.

        Assuming the order instead of reading it swaps two channels, and the
        trace then disagrees with its own base calls at every position.
        """
        for order in ("GATC", "ACGT", "TACG"):
            trace = read_ab1(build_ab1(READ, order=order))
            # The peak under each called base must be that base's channel.
            for i, base in enumerate(READ[:20]):
                at = trace.peaks[i]
                tallest = max("ACGT", key=lambda b: trace.traces[b][at])
                assert tallest == base, f"{order}: base {i} reads as {tallest}"

    def test_short_values_stored_inside_the_entry_are_read(self):
        """Four bytes or fewer live in the offset field, not at it.

        A reader that always follows the offset gets a pointer where the data
        should be — and FWO_, which decides channel order, is exactly 4 bytes.
        """
        trace = read_ab1(build_ab1("ACGT", [30, 30, 30, 30]))
        assert trace.sequence == "ACGT"

    def test_peaks_locate_each_base_in_the_samples(self):
        trace = read_ab1(build_ab1(READ, spacing=14))
        assert len(trace.peaks) == len(READ)
        assert trace.peaks == sorted(trace.peaks), "peaks must run left to right"


class TestQualityTrimming:
    """Sanger reads are unreliable at both ends, by an amount that varies."""

    def test_noisy_ends_are_cut_and_the_good_middle_kept(self):
        quality = [4] * 20 + [50] * 60 + [4] * 25
        trace = read_ab1(build_ab1("A" * 105, quality))
        start, stop = trace.trim()
        assert 15 <= start <= 21
        assert 79 <= stop <= 85

    def test_a_clean_read_is_kept_whole(self):
        """A fixed trim throws away good sequence on a good trace."""
        trace = read_ab1(build_ab1(READ, [55] * len(READ)))
        assert trace.trim() == (0, len(READ))

    def test_an_unusable_trace_trims_to_nothing(self):
        """Better to report no usable sequence than to compare rubbish."""
        trace = read_ab1(build_ab1("A" * 60, [3] * 60))
        start, stop = trace.trim()
        assert stop - start == 0

    def test_the_cutoff_is_adjustable(self):
        quality = [12] * 40
        trace = read_ab1(build_ab1("A" * 40, quality))
        assert trace.trim(cutoff=DEFAULT_TRIM_QUALITY) == (0, 0)
        assert trace.trim(cutoff=5)[1] - trace.trim(cutoff=5)[0] == 40

    def test_trimmed_sequence_matches_the_trimmed_range(self):
        quality = [4] * 10 + [50] * 40 + [4] * 11
        trace = read_ab1(build_ab1(READ, quality))
        start, stop = trace.trim()
        assert trace.trimmed_sequence() == trace.sequence[start:stop]


class TestSummary:
    def test_reports_what_decides_whether_to_use_a_read(self):
        quality = [5] * 15 + [45] * 40 + [5] * 6
        s = summarise(read_ab1(build_ab1(READ, quality), ))
        assert s.length == len(READ)
        assert s.trimmed_length == s.trim_stop - s.trim_start
        assert s.high_quality_bases <= s.trimmed_length
        assert s.usable is False or s.trimmed_length >= 50

    def test_a_short_read_is_not_usable(self):
        s = summarise(read_ab1(build_ab1("ACGT" * 8, [45] * 32)))
        assert s.usable is False, "32 bases cannot verify a construct"


class TestTraceWindow:
    """The drawing sends samples around one base, not the whole trace."""

    def test_window_covers_the_neighbouring_bases(self):
        trace = read_ab1(build_ab1(READ))
        w = trace.window(20, span=4)
        assert w["centre"] == 20
        assert [b["index"] for b in w["bases"]] == list(range(16, 25))
        assert w["bases"][4]["base"] == READ[20]

    def test_window_positions_are_relative_to_the_slice(self):
        """The client draws from zero; absolute sample numbers would be off
        the left edge of every window but the first."""
        trace = read_ab1(build_ab1(READ))
        w = trace.window(30, span=3)
        first, last = w["samples"]
        assert all(0 <= b["at"] <= last - first for b in w["bases"])
        assert len(next(iter(w["traces"].values()))) == last - first

    def test_window_at_either_end_does_not_run_off(self):
        trace = read_ab1(build_ab1(READ))
        for index in (0, 1, len(READ) - 1):
            w = trace.window(index, span=6)
            assert w["bases"], f"no bases returned at {index}"
            assert all(0 <= b["index"] < len(READ) for b in w["bases"])

    def test_an_out_of_range_index_returns_an_empty_window(self):
        trace = read_ab1(build_ab1(READ))
        assert trace.window(10_000)["bases"] == []


class TestBadFiles:
    """These arrive by email from a facility, and are routinely the wrong file."""

    def test_a_pdf_is_named_as_such(self):
        with pytest.raises(SequenceError, match="not an .ab1"):
            read_ab1(b"%PDF-1.4" + b"\x00" * 400)

    def test_a_truncated_download_is_reported(self):
        with pytest.raises(SequenceError, match="too small"):
            read_ab1(b"ABIF" + b"\x00" * 20)

    def test_a_trace_without_base_calls_says_so(self):
        with pytest.raises(SequenceError, match="no base calls"):
            read_ab1(build_ab1(READ, omit={"PBAS"}))

    def test_a_file_pointing_past_its_own_end_does_not_crash(self):
        """A half-finished download has a valid header and no data."""
        blob = bytearray(build_ab1(READ))
        del blob[len(blob) // 2 :]
        with pytest.raises(SequenceError):
            read_ab1(bytes(blob))

    def test_missing_quality_is_survivable(self):
        """Some older instruments write no PCON. The bases are still usable;
        only the confidence is unknown, and trimming must not then guess."""
        trace = read_ab1(build_ab1(READ, omit={"PCON"}))
        assert trace.sequence == READ
        assert trace.quality == []
        assert trace.trim() == (0, len(READ))


class TestVerifyingAgainstTheTrace:
    """The reason for reading the trace at all.

    A read differs from the design at one base. Whether that is a mutation
    worth reordering an oligo over, or the basecaller picking between a peak
    and its neighbour's shoulder, is not in the letters. It is in the peak.
    """

    DESIGN = (
        "ATGGCTAGCAAAGAACTGGTTACCGCTCTGTATCTGGTGTGCGGCGAACGCGGCTTTTTCTACACCCCG"
        "AAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGGGCGGCGGCCCGGGCGCGGGC"
        "AGCCTGCAGCCGCTGGCGCTGGAAGGCAGCCTGCAGAAACGCGGCATCGTGGAACAGTGCTGCACCAGC"
        "ATCTGCAGCCTGTACCAGCTGGAAAACTACTGCAACGGCGGCTTTGTGAACCAGCATCTGTGCGGCAGC"
    )

    def _read_with_a_change(self, at: int = 100):
        """A 200-base read off the design, with one substituted base."""
        from gsynth_engine.verify import verify_read

        piece = list(self.DESIGN[30:230])
        piece[at] = {"A": "G", "G": "A", "C": "T", "T": "C"}[piece[at]]
        return "".join(piece), verify_read

    def test_a_change_on_a_clean_peak_is_reported_as_real(self):
        read, verify_read = self._read_with_a_change()
        trace = read_ab1(build_ab1(read, [48] * len(read)))

        result = verify_read(self.DESIGN, read, trace=trace)
        assert len(result.differences) == 1
        difference = result.differences[0]
        assert difference.quality == 48
        assert difference.confident is True
        assert "check the trace" not in difference.description

    def test_the_same_change_on_a_poor_peak_is_flagged_instead(self):
        """Identical letters, opposite conclusion — which is the whole point."""
        read, verify_read = self._read_with_a_change()
        quality = [48] * len(read)
        quality[100] = 8
        trace = read_ab1(build_ab1(read, quality))

        result = verify_read(self.DESIGN, read, trace=trace)
        difference = result.differences[0]
        assert difference.quality == 8
        assert difference.confident is False
        assert "Q8" in difference.description
        assert result.unconfident_differences == [difference]
        assert result.confirmed_differences == []

    def test_a_read_given_as_letters_claims_no_confidence(self):
        """Without a trace there is nothing to judge, and guessing would be
        worse than saying so — an unmarked difference must mean 'unknown',
        not 'fine'."""
        read, verify_read = self._read_with_a_change()
        result = verify_read(self.DESIGN, read)
        difference = result.differences[0]
        assert difference.quality is None
        assert difference.confident is None
        assert result.confirmed_differences == [difference]

    def test_a_reverse_read_scores_the_base_it_actually_called(self):
        """The trace is in the read's own direction; the comparison is not.

        Mapping the difference back without undoing the flip reads the
        quality of the base the same distance from the *other* end — which
        is a real number, from the wrong base, and looks entirely plausible.
        """
        from gsynth_engine.sequence import reverse_complement
        from gsynth_engine.verify import verify_read

        forward, _ = self._read_with_a_change(at=100)
        read = reverse_complement(forward)
        bad_at = len(read) - 1 - 100            # the same base, counted the other way

        quality = [48] * len(read)
        quality[bad_at] = 6
        trace = read_ab1(build_ab1(read, quality))

        result = verify_read(self.DESIGN, read, trace=trace)
        assert result.reverse_complemented
        assert len(result.differences) == 1
        assert result.differences[0].read_index == bad_at
        assert result.differences[0].quality == 6

    def test_a_trace_trims_by_quality_not_by_a_fixed_count(self):
        """A clean read keeps its ends; the fixed 30 would discard 60 good
        bases, and on a bad read it would keep 30 useless ones."""
        from gsynth_engine.verify import verify_read

        read = self.DESIGN[30:230]
        trace = read_ab1(build_ab1(read, [50] * len(read)))

        by_quality = verify_read(self.DESIGN, read, trace=trace)
        by_count = verify_read(self.DESIGN, read, trim=30)
        assert by_quality.trimmed_start == 0
        assert by_count.trimmed_start == 30
        assert by_quality.covered > by_count.covered

    def test_the_report_carries_the_traces_mean_quality(self):
        from gsynth_engine.verify import verify

        read = self.DESIGN[30:230]
        trace = read_ab1(build_ab1(read, [42] * len(read)), name="fwd")
        report = verify(self.DESIGN, {"fwd": read}, traces={"fwd": trace})
        assert report.reads[0].mean_quality == 42.0
