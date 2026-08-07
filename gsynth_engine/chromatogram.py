"""Sanger trace files: the base calls, and how much to believe each one.

The facility returns an `.ab1`. Until now the workflow was to open it in
another program, read the bases off, and paste them in — which throws away
the only part that answers the question people actually have. A read differs
from the design at one position: is that a mutation, or a bad call?

The trace answers it and the letters do not. A clean base is one sharp peak
at Q50 or better; an artefact is a shoulder on its neighbour at Q10, and both
arrive as the same letter. So this module keeps the quality values and the
four channels, and `verify.py` uses them to say which differences are worth a
second look and which are noise.

ABIF is a documented format — a header, a directory of tagged entries, and
their data — so it is read here rather than pulled in with a library. The
engine has no runtime dependencies and this is not the place to acquire one.
The spec is Applied Biosystems' "ABIF File Format" (2006, rev. 2009).
"""
from __future__ import annotations

import struct
from dataclasses import dataclass

from gsynth_engine.sequence import SequenceError

#: Everything past this is a directory of tags; the header is fixed.
_MAGIC = b"ABIF"
# name, number, type, element size, count, byte count, offset, handle.
# The offset field is kept as raw bytes: when the data is four bytes or
# fewer it is stored *in* that field rather than pointed to by it.
_ENTRY = ">4sihhii4si"
_ENTRY_LEN = 28

#: struct codes for the element types this reader needs. The format defines
#: more (dates, times, user-defined), none of which carry trace data.
_TYPES: dict[int, tuple[str, int]] = {
    1: ("B", 1),      # byte
    2: ("c", 1),      # char
    3: ("H", 2),      # word
    4: ("h", 2),      # short
    5: ("i", 4),      # long
    7: ("f", 4),      # float
    8: ("d", 8),      # double
}

#: A trace longer than this is not a Sanger read; it is a mistake or an
#: attack. A 1.2 kb read is about 15 000 samples.
MAX_SAMPLES = 200_000

#: Mott trimming's default: bases worse than Q13 are more likely wrong than
#: a coin toss over a 20-base window. Phred's own default, and Sequencher's.
DEFAULT_TRIM_QUALITY = 13


@dataclass(frozen=True)
class Chromatogram:
    """One Sanger read: what was called, and how well."""

    sequence: str                       #: base calls, 5'→3'
    quality: list[int]                  #: Phred score per called base
    peaks: list[int]                    #: sample index of each base's peak
    traces: dict[str, list[int]]        #: "A"/"C"/"G"/"T" → signal per sample
    name: str = ""

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sample_count(self) -> int:
        return len(next(iter(self.traces.values()), []))

    @property
    def mean_quality(self) -> float:
        return sum(self.quality) / len(self.quality) if self.quality else 0.0

    def quality_at(self, index: int) -> int:
        """The Phred score of one called base, or 0 past either end."""
        return self.quality[index] if 0 <= index < len(self.quality) else 0

    def trim(self, cutoff: int = DEFAULT_TRIM_QUALITY) -> tuple[int, int]:
        """The good stretch, by Richard Mott's algorithm.

        Returns a half-open ``(start, stop)`` over the called bases.

        Sanger reads are unreliable at both ends — primer-proximal noise for
        the first 20-40 bases, decaying signal after 700-900. Cutting a fixed
        number off each end, as this workflow did before, is wrong in both
        directions on the same read: it keeps rubbish on a bad trace and
        discards good sequence on a clean one.

        Each base scores ``cutoff_probability - its own error probability``,
        so good bases score positive and bad ones negative; the best-scoring
        run is the region to keep.
        """
        if not self.quality:
            return 0, len(self.sequence)

        limit = 10 ** (-cutoff / 10)
        best = best_start = best_stop = 0
        run = run_start = 0

        for i, q in enumerate(self.quality):
            run += limit - 10 ** (-q / 10)
            if run < 0:                       # this run cannot start a better one
                run, run_start = 0.0, i + 1
            elif run > best:
                best, best_start, best_stop = run, run_start, i + 1

        # An entirely poor trace has no positive-scoring run at all.
        return (best_start, best_stop) if best > 0 else (0, 0)

    def trimmed_sequence(self, cutoff: int = DEFAULT_TRIM_QUALITY) -> str:
        start, stop = self.trim(cutoff)
        return self.sequence[start:stop]

    def window(self, index: int, span: int = 5) -> dict:
        """The trace around one called base, for drawing.

        Sending every sample of every channel to a browser is megabytes for
        a picture of ten bases. This returns just the samples spanned by the
        neighbouring peaks, with the base positions marked inside them.
        """
        if not self.peaks or not 0 <= index < len(self.peaks):
            return {"samples": [], "traces": {}, "bases": [], "centre": 0}

        low = max(0, index - span)
        high = min(len(self.peaks) - 1, index + span)

        # Half a peak-spacing of margin, so the outer peaks are not clipped.
        spacing = self._spacing()
        first = max(0, self.peaks[low] - spacing // 2)
        last = min(self.sample_count, self.peaks[high] + spacing // 2 + 1)

        return {
            "samples": [first, last],
            "traces": {b: t[first:last] for b, t in self.traces.items()},
            "bases": [
                {
                    "index": i,
                    "base": self.sequence[i],
                    "quality": self.quality_at(i),
                    "at": self.peaks[i] - first,
                }
                for i in range(low, high + 1)
                if i < len(self.sequence)
            ],
            "centre": index,
        }

    def _spacing(self) -> int:
        """Typical samples between adjacent peaks — about 12 on a 3730."""
        if len(self.peaks) < 2:
            return 12
        span = self.peaks[-1] - self.peaks[0]
        return max(1, span // (len(self.peaks) - 1))


@dataclass(frozen=True)
class _Entry:
    tag: str
    number: int
    kind: int
    count: int
    nbytes: int
    offset: bytes          #: raw; an offset, or the data itself when ≤ 4 bytes


def _read_entry(blob: bytes, at: int) -> _Entry:
    name, number, kind, _size, count, nbytes, offset, _handle = struct.unpack(
        _ENTRY, blob[at : at + _ENTRY_LEN]
    )
    return _Entry(name.decode("ascii", "replace"), number, kind, count, nbytes, offset)


def _value(blob: bytes, entry: _Entry):
    """One entry's data, decoded. Four bytes or fewer are stored in place."""
    kind, count, nbytes = entry.kind, entry.count, entry.nbytes
    if nbytes <= 4:
        raw = entry.offset[:nbytes]
    else:
        start = struct.unpack(">i", entry.offset)[0]
        raw = blob[start : start + nbytes]
    if len(raw) < nbytes:
        raise SequenceError("This trace file is truncated — the data it "
                            "points to is not there. Re-export it from the "
                            "sequencing facility's own file.")
    if kind == 2:                                   # char array: bases, FWO_
        return raw.decode("ascii", "replace")
    if kind in (18, 19):                            # pString / cString
        return raw[1:].decode("ascii", "replace") if kind == 18 else \
               raw.rstrip(b"\x00").decode("ascii", "replace")
    code, width = _TYPES.get(kind, ("", 0))
    if not code:
        return raw
    return list(struct.unpack(f">{count}{code}", raw[: count * width]))


def read_ab1(data: bytes, *, name: str = "") -> Chromatogram:
    """Parse an ABIF trace file.

    Raises `SequenceError` with something the user can act on — these files
    arrive by email from a facility and are routinely the wrong file, a
    renamed `.pdf`, or half a download.
    """
    if len(data) < 128:
        raise SequenceError(
            "That file is too small to be a trace. An .ab1 from a sequencing "
            "facility is normally 100-500 kB."
        )
    if data[:4] != _MAGIC:
        raise SequenceError(
            "That is not an .ab1 trace file — it does not start with the ABIF "
            "marker. If the facility sent a .zip or a .pdf, send the .ab1 "
            "inside it instead."
        )

    try:
        header = _read_entry(data, 6)
        directory = struct.unpack(">i", header.offset)[0]
    except struct.error as exc:                         # pragma: no cover
        raise SequenceError("This trace file's header is damaged.") from exc

    tags: dict[tuple[str, int], object] = {}
    for i in range(header.count):
        at = directory + i * _ENTRY_LEN
        if at + _ENTRY_LEN > len(data):
            break
        entry = _read_entry(data, at)
        try:
            tags[(entry.tag, entry.number)] = _value(data, entry)
        except (struct.error, SequenceError):
            continue                                    # a tag we cannot use

    def first(tag: str, *numbers: int):
        for n in numbers:
            if (tag, n) in tags:
                return tags[(tag, n)]
        return None

    # PBAS2 is the basecaller's own sequence; PBAS1 may have been edited by
    # hand in a viewer. Prefer the edited one when it exists — someone looked.
    sequence = first("PBAS", 1, 2)
    if not isinstance(sequence, str) or not sequence:
        raise SequenceError(
            "This trace has no base calls in it. The facility may have sent "
            "the raw instrument file rather than the analysed one."
        )
    sequence = sequence.upper().replace("-", "N")

    quality = first("PCON", 1, 2) or []
    if isinstance(quality, (bytes, str)):
        quality = list(bytes(quality, "latin-1") if isinstance(quality, str) else quality)
    quality = [int(q) for q in quality][: len(sequence)]

    peaks = [int(p) for p in (first("PLOC", 1, 2) or [])][: len(sequence)]

    # DATA9-12 are the processed channels; FWO_ says which base each carries.
    order = first("FWO_", 1) or "GATC"
    order = "".join(b for b in str(order).upper() if b in "ACGT") or "GATC"
    traces: dict[str, list[int]] = {}
    for i, base in enumerate(order[:4]):
        channel = first("DATA", 9 + i)
        if isinstance(channel, list):
            if len(channel) > MAX_SAMPLES:
                raise SequenceError(
                    f"This trace has {len(channel):,} samples, which is far "
                    "longer than a Sanger read. Check it is a single read and "
                    "not a concatenated file."
                )
            traces[base] = [int(v) for v in channel]

    for base in "ACGT":
        traces.setdefault(base, [])

    # A half-finished download still parses: the base calls sit near the front
    # of the file and the channels at the back. Letting that through hands
    # back letters with no peaks under them, which is the one thing this
    # module exists to prevent — and it would look like a working read.
    carrying = sorted(b for b, t in traces.items() if t)
    if len(carrying) < 4:
        missing = ", ".join(b for b in "ACGT" if b not in carrying)
        raise SequenceError(
            f"This trace file is missing the {missing} channel"
            f"{'s' if len(carrying) < 3 else ''} — a Sanger read has all four. "
            "That means a partial download or the raw instrument file; ask "
            "the facility to send the analysed .ab1 again."
        )
    shortest = min(len(traces[b]) for b in carrying)
    if peaks and max(peaks) >= shortest:
        raise SequenceError(
            f"This trace file stops before its last base: the peaks run to "
            f"sample {max(peaks):,} but the signal ends at {shortest:,}. The "
            "download is incomplete — ask the facility to send it again."
        )

    return Chromatogram(
        sequence=sequence,
        quality=quality,
        peaks=peaks,
        traces=traces,
        name=name,
    )


@dataclass(frozen=True)
class TraceSummary:
    """What a person needs to decide whether to trust a read."""

    name: str
    length: int
    mean_quality: float
    trim_start: int
    trim_stop: int
    trimmed_length: int
    high_quality_bases: int          #: Q20 or better, within the trimmed region
    sample_count: int

    @property
    def usable(self) -> bool:
        """Enough good sequence to be worth comparing to a design at all."""
        return self.trimmed_length >= 50 and self.mean_quality >= 15


def summarise(trace: Chromatogram,
              cutoff: int = DEFAULT_TRIM_QUALITY) -> TraceSummary:
    start, stop = trace.trim(cutoff)
    kept = trace.quality[start:stop]
    return TraceSummary(
        name=trace.name,
        length=trace.length,
        mean_quality=round(trace.mean_quality, 1),
        trim_start=start,
        trim_stop=stop,
        trimmed_length=stop - start,
        high_quality_bases=sum(1 for q in kept if q >= 20),
        sample_count=trace.sample_count,
    )
