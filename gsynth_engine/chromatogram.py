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

ABIF and SCF are documented formats, so both are read here rather than pulled
in with a library.  Some sequencing software exports SCF data without changing
the original ``.ab1`` filename, which is why format detection uses the file's
magic bytes instead of its extension.  The engine has no runtime dependencies
and this is not the place to acquire one.  The specifications are Applied
Biosystems' "ABIF File Format" (2006, rev. 2009) and Staden's SCF v2/v3 format.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass

from gsynth_engine.sequence import SequenceError, reverse_complement

#: Everything past this is a directory of tags; the header is fixed.
_MAGIC = b"ABIF"
_SCF_MAGIC = b".scf"
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

# A Sanger trace normally contains fewer than 2,000 calls.  This deliberately
# generous ceiling also bounds allocations before any section is unpacked.
MAX_BASES = 100_000

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

    def alignment_track(self, start: int, stop: int, *, reverse: bool = False) -> dict:
        """The quality-trimmed evidence needed for a reference-aligned viewer.

        Unlike :meth:`window`, this covers the complete part of a read that
        was actually admitted to verification.  It deliberately excludes the
        discarded noisy ends: showing those bases on the reference track
        would visually imply evidence that the algorithm did not use.

        Reverse reads are returned in reference orientation.  Reversing a
        chromatogram means both reversing sample order *and* complementing
        its channels (A<->T and C<->G); doing only one of those makes the
        coloured peaks disagree with the displayed base calls.
        """
        start = max(0, min(start, self.length))
        stop = max(start, min(stop, self.length))
        if start == stop or not self.peaks:
            return {
                "sequence": "", "qualities": [], "peaks": [],
                "sample_count": 0, "traces": {base: [] for base in "ACGT"},
            }

        spacing = self._spacing()
        first = max(0, self.peaks[start] - spacing // 2)
        last = min(self.sample_count, self.peaks[stop - 1] + spacing // 2 + 1)
        sequence = self.sequence[start:stop]
        qualities = self.quality[start:stop]
        peaks = [peak - first for peak in self.peaks[start:stop]]
        channels = {base: values[first:last] for base, values in self.traces.items()}

        if reverse:
            complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
            width = last - first
            sequence = reverse_complement(sequence)
            qualities = list(reversed(qualities))
            peaks = [width - 1 - peak for peak in reversed(peaks)]
            channels = {
                base: list(reversed(channels.get(complement[base], [])))
                for base in "ACGT"
            }

        return {
            "sequence": sequence,
            "qualities": qualities,
            "peaks": peaks,
            "sample_count": last - first,
            "traces": channels,
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


def _section(data: bytes, offset: int, size: int, description: str) -> bytes:
    """Return a bounded SCF section, rejecting wraparound and truncation."""
    if offset < 128 or size < 0 or offset > len(data) or size > len(data) - offset:
        raise SequenceError(
            f"This SCF trace is truncated or damaged: its {description} section "
            "points outside the file. Re-export the analysed trace from the "
            "sequencing software."
        )
    return data[offset : offset + size]


def _delta_delta_decode(values: list[int], modulus: int) -> list[int]:
    """Undo SCF v3's two successive modulo-delta encodings."""
    decoded = list(values)
    for _ in range(2):
        running = 0
        for index, value in enumerate(decoded):
            running = (running + value) % modulus
            decoded[index] = running
    return decoded


def read_scf(data: bytes, *, name: str = "") -> Chromatogram:
    """Parse a Staden SCF v2 or v3 chromatogram.

    SCF stores the same scientific evidence as ABIF: four electropherogram
    channels, called bases, peak locations, and per-base probabilities.  SCF
    v3 groups channels and base fields into blocks and delta-delta encodes the
    samples; v2 stores interleaved sample and base records.
    """
    if len(data) < 128:
        raise SequenceError(
            "That file is too small to be a Sanger trace. An ABIF or SCF file "
            "from a sequencing facility is normally several kilobytes."
        )
    if data[:4] != _SCF_MAGIC:
        raise SequenceError(
            "That is not an SCF trace file — it does not start with the .scf "
            "marker. Upload the analysed ABIF or SCF chromatogram, not a PDF "
            "or sequence-only export."
        )

    try:
        (
            _magic, samples, samples_offset, bases, _left_clip, _right_clip,
            bases_offset, comments_size, comments_offset, version_raw,
            sample_size, _code_set, private_size, private_offset, *_spare,
        ) = struct.unpack(">4s8I4s4I18I", data[:128])
    except struct.error as exc:  # pragma: no cover - length checked above
        raise SequenceError("This SCF trace file's header is damaged.") from exc

    try:
        version_text = version_raw.decode("ascii")
        major = int(version_text.split(".", 1)[0])
    except (UnicodeDecodeError, ValueError) as exc:
        raise SequenceError("This SCF trace has an invalid format version.") from exc
    if major not in (2, 3):
        raise SequenceError(
            f"SCF version {version_text!r} is not supported. Re-export the "
            "trace as SCF v3 or ABIF."
        )
    if samples == 0 or samples > MAX_SAMPLES:
        raise SequenceError(
            f"This trace declares {samples:,} samples, outside the supported "
            f"range of 1–{MAX_SAMPLES:,}. Check that it is one Sanger read."
        )
    if bases == 0 or bases > MAX_BASES:
        raise SequenceError(
            f"This trace declares {bases:,} base calls, outside the supported "
            f"range of 1–{MAX_BASES:,}. Check that it is one Sanger read."
        )
    if sample_size not in (1, 2):
        raise SequenceError(
            f"This SCF trace uses {sample_size}-byte samples; the format allows "
            "only one- or two-byte samples."
        )

    sample_bytes = samples * 4 * sample_size
    raw_samples = _section(data, samples_offset, sample_bytes, "sample")
    code = "B" if sample_size == 1 else "H"
    width = sample_size
    traces: dict[str, list[int]] = {base: [] for base in "ACGT"}

    if major == 3:
        for channel_index, base in enumerate("ACGT"):
            start = channel_index * samples * width
            raw = raw_samples[start : start + samples * width]
            encoded = list(struct.unpack(f">{samples}{code}", raw))
            traces[base] = _delta_delta_decode(encoded, 1 << (8 * width))
    else:
        values = struct.unpack(f">{samples * 4}{code}", raw_samples)
        for channel_index, base in enumerate("ACGT"):
            traces[base] = [
                int(values[sample_index * 4 + channel_index])
                for sample_index in range(samples)
            ]

    base_bytes = bases * 12
    raw_bases = _section(data, bases_offset, base_bytes, "base-call")

    # Validate optional sections too.  They are not needed for verification,
    # but a bogus size must not let a damaged file masquerade as a valid one.
    if comments_size:
        _section(data, comments_offset, comments_size, "comments")
    if private_size:
        _section(data, private_offset, private_size, "private data")
    peaks: list[int] = []
    probabilities: dict[str, list[int]] = {base: [] for base in "ACGT"}

    if major == 3:
        peaks = list(struct.unpack(f">{bases}I", raw_bases[: bases * 4]))
        cursor = bases * 4
        for base in "ACGT":
            probabilities[base] = list(raw_bases[cursor : cursor + bases])
            cursor += bases
        calls = raw_bases[cursor : cursor + bases]
    else:
        calls_buffer = bytearray()
        for index in range(bases):
            at = index * 12
            peak, pa, pc, pg, pt, call, _sub, _ins, _del = struct.unpack(
                ">I4B1s3B", raw_bases[at : at + 12]
            )
            peaks.append(peak)
            for base, probability in zip("ACGT", (pa, pc, pg, pt), strict=True):
                probabilities[base].append(probability)
            calls_buffer += call
        calls = bytes(calls_buffer)

    sequence = calls.decode("ascii", "replace").upper().replace("-", "N")
    sequence = "".join(base if base in "ACGTN" else "N" for base in sequence)
    quality = [
        probabilities[call][index]
        if call in probabilities
        else max(probabilities[base][index] for base in "ACGT")
        for index, call in enumerate(sequence)
    ]

    if len(peaks) != len(sequence) or any(
        peak >= samples or (index and peak < peaks[index - 1])
        for index, peak in enumerate(peaks)
    ):
        raise SequenceError(
            "This SCF trace has invalid or out-of-order peak locations. "
            "Re-export the analysed trace from the sequencing software."
        )

    return Chromatogram(
        sequence=sequence,
        quality=quality,
        peaks=[int(peak) for peak in peaks],
        traces=traces,
        name=name,
    )


def read_trace(data: bytes, *, name: str = "") -> Chromatogram:
    """Read a Sanger chromatogram by its contents, not its filename."""
    magic = data[:4]
    if magic == _MAGIC:
        return read_ab1(data, name=name)
    if magic == _SCF_MAGIC:
        return read_scf(data, name=name)
    if len(data) < 128:
        raise SequenceError(
            "That file is too small to be a Sanger trace. Upload the analysed "
            "ABIF (.ab1) or SCF chromatogram from the sequencing facility."
        )
    raise SequenceError(
        "That is not a supported Sanger trace: it has neither the ABIF nor "
        "the SCF marker. Upload the analysed .ab1 or .scf file, not a PDF, "
        "archive, or sequence-only export."
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
