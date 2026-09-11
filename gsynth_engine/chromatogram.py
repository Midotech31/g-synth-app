from __future__ import annotations

import struct
from dataclasses import dataclass

from gsynth_engine.sequence import SequenceError, reverse_complement

_MAGIC = b"ABIF"
_SCF_MAGIC = b".scf"


_ENTRY = ">4sihhii4si"
_ENTRY_LEN = 28


_TYPES: dict[int, tuple[str, int]] = {
    1: ("B", 1),
    2: ("c", 1),
    3: ("H", 2),
    4: ("h", 2),
    5: ("i", 4),
    7: ("f", 4),
    8: ("d", 8),
}


MAX_SAMPLES = 200_000


MAX_BASES = 100_000


DEFAULT_TRIM_QUALITY = 13


@dataclass(frozen=True)
class Chromatogram:


    sequence: str
    quality: list[int]
    peaks: list[int]
    traces: dict[str, list[int]]
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

        return self.quality[index] if 0 <= index < len(self.quality) else 0

    def trim(self, cutoff: int = DEFAULT_TRIM_QUALITY) -> tuple[int, int]:

        if not self.quality:
            return 0, len(self.sequence)

        limit = 10 ** (-cutoff / 10)
        best = best_start = best_stop = 0
        run = run_start = 0

        for i, q in enumerate(self.quality):
            run += limit - 10 ** (-q / 10)
            if run < 0:
                run, run_start = 0.0, i + 1
            elif run > best:
                best, best_start, best_stop = run, run_start, i + 1


        return (best_start, best_stop) if best > 0 else (0, 0)

    def trimmed_sequence(self, cutoff: int = DEFAULT_TRIM_QUALITY) -> str:
        start, stop = self.trim(cutoff)
        return self.sequence[start:stop]

    def window(self, index: int, span: int = 5) -> dict:

        if not self.peaks or not 0 <= index < len(self.peaks):
            return {"samples": [], "traces": {}, "bases": [], "centre": 0}

        low = max(0, index - span)
        high = min(len(self.peaks) - 1, index + span)


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
    offset: bytes


def _read_entry(blob: bytes, at: int) -> _Entry:
    name, number, kind, _size, count, nbytes, offset, _handle = struct.unpack(
        _ENTRY, blob[at : at + _ENTRY_LEN]
    )
    return _Entry(name.decode("ascii", "replace"), number, kind, count, nbytes, offset)


def _value(blob: bytes, entry: _Entry):

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
    if kind == 2:
        return raw.decode("ascii", "replace")
    if kind in (18, 19):
        return raw[1:].decode("ascii", "replace") if kind == 18 else \
               raw.rstrip(b"\x00").decode("ascii", "replace")
    code, width = _TYPES.get(kind, ("", 0))
    if not code:
        return raw
    return list(struct.unpack(f">{count}{code}", raw[: count * width]))


def read_ab1(data: bytes, *, name: str = "") -> Chromatogram:

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

    header = _read_entry(data, 6)
    directory = struct.unpack(">i", header.offset)[0]

    tags: dict[tuple[str, int], object] = {}
    for i in range(header.count):
        at = directory + i * _ENTRY_LEN
        if at + _ENTRY_LEN > len(data):
            break
        entry = _read_entry(data, at)
        try:
            tags[(entry.tag, entry.number)] = _value(data, entry)
        except (struct.error, SequenceError):
            continue

    def first(tag: str, *numbers: int):
        for n in numbers:
            if (tag, n) in tags:
                return tags[(tag, n)]
        return None


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

    if offset < 128 or size < 0 or offset > len(data) or size > len(data) - offset:
        raise SequenceError(
            f"This SCF trace is truncated or damaged: its {description} section "
            "points outside the file. Re-export the analysed trace from the "
            "sequencing software."
        )
    return data[offset : offset + size]


def _delta_delta_decode(values: list[int], modulus: int) -> list[int]:

    decoded = list(values)
    for _ in range(2):
        running = 0
        for index, value in enumerate(decoded):
            running = (running + value) % modulus
            decoded[index] = running
    return decoded


def read_scf(data: bytes, *, name: str = "") -> Chromatogram:

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

    (
        _magic, samples, samples_offset, bases, _left_clip, _right_clip,
        bases_offset, comments_size, comments_offset, version_raw,
        sample_size, _code_set, private_size, private_offset, *_spare,
    ) = struct.unpack(">4s8I4s4I18I", data[:128])

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


    name: str
    length: int
    mean_quality: float
    trim_start: int
    trim_stop: int
    trimmed_length: int
    high_quality_bases: int
    sample_count: int

    @property
    def usable(self) -> bool:

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
