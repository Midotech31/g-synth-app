from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from gsynth_engine.chromatogram import DEFAULT_TRIM_QUALITY, Chromatogram
from gsynth_engine.cloning import translate
from gsynth_engine.sequence import SequenceError, clean_dna, reverse_complement

ANCHOR = 14


BAND = 24


@dataclass(frozen=True)
class Difference:


    kind: str
    position: int
    expected: str
    found: str

    codon: int | None = None
    residue: int | None = None
    from_residue: str = ""
    to_residue: str = ""
    silent: bool | None = None


    read_index: int | None = None
    quality: int | None = None

    @property
    def confident(self) -> bool | None:

        return None if self.quality is None else self.quality >= 20

    @property
    def description(self) -> str:
        where = f"position {self.position + 1}"
        if self.kind == "substitution":
            change = f"{self.expected}→{self.found} at {where}"
        elif self.kind == "insertion":
            change = f"{self.found} inserted at {where}"
        else:
            change = f"{self.expected} missing at {where}"
        if self.quality is not None and self.quality < 20:
            change += f" (Q{self.quality} — check the trace)"
        if self.residue is None:
            return change
        if self.silent:
            return f"{change} — silent, residue {self.residue} stays {self.from_residue}"
        return (
            f"{change} — residue {self.residue} "
            f"{self.from_residue}→{self.to_residue or '?'}"
        )


@dataclass
class ReadAlignment:


    name: str
    length: int

    start: int
    end: int
    reverse_complemented: bool
    identity: float
    matched: int
    differences: list[Difference] = field(default_factory=list)
    trimmed_start: int = 0
    trimmed_end: int = 0
    warnings: list[str] = field(default_factory=list)

    mean_quality: float | None = None

    @property
    def unconfident_differences(self) -> list[Difference]:

        return [d for d in self.differences if d.confident is False]

    @property
    def confirmed_differences(self) -> list[Difference]:

        return [d for d in self.differences if d.confident is not False]

    @property
    def covered(self) -> int:

        return self.end - self.start

    @property
    def is_clean(self) -> bool:

        return not self.differences


@dataclass
class VerificationReport:


    design_length: int
    reads: list[ReadAlignment] = field(default_factory=list)

    gaps: list[tuple[int, int]] = field(default_factory=list)
    coverage: float = 0.0
    warnings: list[str] = field(default_factory=list)

    @property
    def differences(self) -> list[Difference]:
        seen: dict[tuple[str, int, str], Difference] = {}
        for read in self.reads:
            for difference in read.differences:
                seen[(difference.kind, difference.position, difference.found)] = difference
        return sorted(seen.values(), key=lambda d: d.position)

    @property
    def is_verified(self) -> bool:

        return bool(self.reads) and self.fully_covered and not self.differences

    @property
    def fully_covered(self) -> bool:

        return not self.gaps


@dataclass(frozen=True)
class ConsensusPosition:


    position: int
    reference: str
    call: str
    supporting_reads: tuple[str, ...]
    qualities: tuple[int, ...]
    combined_quality: int
    agreement: bool


@dataclass
class ConsensusReport:


    reference_length: int
    sequence: str
    positions: list[ConsensusPosition] = field(default_factory=list)
    gaps: list[tuple[int, int]] = field(default_factory=list)
    coverage: float = 0.0
    identity: float = 0.0
    bidirectional_overlap: float = 0.0
    bidirectional_agreement: float = 0.0
    warnings: list[str] = field(default_factory=list)

    @property
    def fully_covered(self) -> bool:
        return not self.gaps

    @property
    def differences(self) -> list[ConsensusPosition]:
        return [
            position for position in self.positions
            if position.call not in {"N", position.reference}
        ]


def _anchors(sequence: str, k: int = ANCHOR) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i in range(len(sequence) - k + 1):
        index.setdefault(sequence[i : i + k], []).append(i)
    return index


def _locate(design: str, read: str, *, circular: bool) -> tuple[int, bool] | None:

    scan = design + design[: ANCHOR - 1] if circular else design
    index = _anchors(scan)

    best: tuple[int, int, bool] | None = None
    for flipped, candidate in ((False, read), (True, reverse_complement(read))):
        votes: Counter[int] = Counter()
        for q in range(0, len(candidate) - ANCHOR + 1, 3):
            positions = index.get(candidate[q : q + ANCHOR])
            if positions and len(positions) == 1:
                votes[positions[0] - q] += 1
        if not votes:
            continue
        offset, count = votes.most_common(1)[0]
        if best is None or count > best[0]:
            best = (count, offset, flipped)

    if best is None or best[0] < 2:
        return None
    return best[1], best[2]


def _align(
    design: str, read: str, offset: int, *, circular: bool,
) -> tuple[
    list[tuple[int, str, str, int]],
    int,
    list[tuple[int | None, str, str, int | None]],
]:

    n = len(read)
    length = len(design)
    if circular:
        window_start = offset - BAND
        m = len(read) + 2 * BAND
        window = "".join(design[(window_start + i) % length] for i in range(m))
    else:
        window_start = max(0, offset - BAND)
        window_stop = min(length, offset + len(read) + BAND)
        window = design[window_start:window_stop]
        m = len(window)

    MATCH, MISMATCH, GAP = 1, -1, -2
    OUTSIDE = -(1 << 30)
    WIDTH = 2 * BAND + 4


    previous = [0] * (m + 1)
    current = [OUTSIDE] * (m + 1)
    trace: list[bytearray] = []
    bands: list[int] = []

    for i in range(1, n + 1):
        low = max(1, i - 1)
        high = min(m, i + 2 * BAND + 1)


        current[0] = i * GAP
        for j in range(max(1, low - 1), min(m, high + 2) + 1):
            current[j] = OUTSIDE

        row = bytearray(WIDTH)
        residue = read[i - 1]

        for j in range(low, high + 1):
            diagonal = previous[j - 1] + (
                MATCH if residue == window[j - 1] else MISMATCH
            )
            up = previous[j] + GAP
            left = current[j - 1] + GAP
            best, step = diagonal, 0
            if up > best:
                best, step = up, 1
            if left > best:
                best, step = left, 2
            current[j] = best
            row[j - low] = step

        trace.append(row)
        bands.append(low)
        previous, current = current, previous


    last_low = bands[-1] if bands else 1
    last_high = min(m, n + 2 * BAND + 1)
    j = max(range(last_low, last_high + 1), key=lambda x: previous[x]) if n else 0
    i = n
    operations: list[tuple[int, str, str, int]] = []
    columns: list[tuple[int | None, str, str, int | None]] = []
    matched = 0

    while i > 0:
        low = bands[i - 1]
        if not low <= j < low + WIDTH:
            raise SequenceError("The read alignment exceeds the supported alignment window; it cannot be verified completely.")
        step = trace[i - 1][j - low]
        if step == 0:
            expected, found = window[j - 1], read[i - 1]
            position = _reference_position(window_start + j - 1, length, circular)
            columns.append((position, expected, found, i - 1))
            if expected == found:
                matched += 1
            else:
                operations.append(
                    (position, expected, found, i - 1))
            i, j = i - 1, j - 1
        elif step == 1:
            columns.append((None, "", read[i - 1], i - 1))
            operations.append((
                _reference_position(window_start + j, length, circular),
                "", read[i - 1], i - 1,
            ))
            i -= 1
        else:
            position = _reference_position(window_start + j - 1, length, circular)
            columns.append((position, window[j - 1], "", None))
            operations.append(
                (position, window[j - 1], "", max(0, i - 1)))
            j -= 1

    operations.reverse()
    columns.reverse()
    return operations, matched, columns


def _reference_position(position: int, length: int, circular: bool) -> int:

    return position % length if circular else min(max(position, 0), length)


def _trim(read: str, *, trim: int) -> tuple[str, int, int]:

    if trim <= 0 or len(read) <= 2 * trim:
        return read, 0, 0
    return read[trim : len(read) - trim], trim, trim


def _coding_effect(
    design: str, position: int, found: str, coding_start: int, coding_end: int,
) -> dict:

    if not (coding_start <= position < coding_end):
        return {}

    offset = position - coding_start
    codon_index = offset // 3
    codon_start = coding_start + codon_index * 3
    if codon_start + 3 > len(design):
        return {}

    original = design[codon_start : codon_start + 3]
    mutated = list(original)
    mutated[offset % 3] = found
    changed = "".join(mutated)

    from_residue = translate(original)
    to_residue = translate(changed)
    return {
        "codon": codon_index + 1,
        "residue": codon_index + 1,
        "from_residue": from_residue,
        "to_residue": to_residue,
        "silent": from_residue == to_residue,
    }


def verify_read(
    design: str,
    read: str,
    *,
    name: str = "read",
    circular: bool = False,
    trim: int = 30,
    coding_start: int | None = None,
    coding_end: int | None = None,
    trace: Chromatogram | None = None,
    trim_quality: int = DEFAULT_TRIM_QUALITY,
) -> ReadAlignment:

    template = clean_dna(design)
    raw = clean_dna(read)
    if not template:
        raise SequenceError("The design is empty.")
    if not raw:
        raise SequenceError(f"{name} is empty.")

    if trace is not None and trace.quality:

        start, stop = trace.trim(trim_quality)
        trimmed = raw[start:stop]
        cut_start, cut_end = start, len(raw) - stop
    else:
        trimmed, cut_start, cut_end = _trim(raw, trim=trim)
    if len(trimmed) < ANCHOR:
        raise SequenceError(
            f"{name} is only {len(raw)} bases; after trimming {trim} from each "
            f"end there is not enough left to place it."
        )

    placed = _locate(template, trimmed, circular=circular)
    if placed is None:
        raise SequenceError(
            f"{name} does not match the design anywhere. Check it is the right "
            f"construct, and that the read is not mostly primer or noise."
        )
    offset, flipped = placed
    full_oriented = reverse_complement(trimmed) if flipped else trimmed


    oriented_prefix = 0
    oriented_suffix = 0
    aligned_offset = offset
    oriented = full_oriented
    if not circular:
        oriented_prefix = max(0, -offset)
        oriented_suffix = max(0, offset + len(full_oriented) - len(template))
        stop = len(full_oriented) - oriented_suffix if oriented_suffix else len(full_oriented)
        oriented = full_oriented[oriented_prefix:stop]
        aligned_offset = max(0, offset)
        if len(oriented) < ANCHOR:
            raise SequenceError(
                f"{name} overlaps only {len(oriented)} bases of the design; "
                "there is not enough reference-aligned sequence to verify it."
            )

    operations, matched, _columns = _align(
        template, oriented, aligned_offset, circular=circular,
    )

    differences: list[Difference] = []
    for position, expected, found, read_at in operations:
        kind = (
            "substitution" if expected and found
            else "insertion" if not expected else "deletion"
        )
        effect = (
            _coding_effect(template, position, found, coding_start, coding_end)
            if kind == "substitution" and coding_start is not None and coding_end is not None
            else {}
        )
        full_oriented_at = oriented_prefix + read_at
        index = (
            len(full_oriented) - 1 - full_oriented_at
            if flipped else full_oriented_at
        )
        index += cut_start
        differences.append(Difference(
            kind=kind, position=position, expected=expected, found=found,
            read_index=index if trace is not None else None,
            quality=trace.quality_at(index) if trace is not None else None,
            **effect,
        ))

    total = matched + len(differences)
    identity = round(100.0 * matched / total, 2) if total else 0.0

    warnings: list[str] = []
    if identity < 98 and differences:
        warnings.append(
            f"{name} matches the design at only {identity:.1f}%. That is more "
            f"disagreement than a good read has — check the trace quality "
            f"before treating these as real changes."
        )

    return ReadAlignment(
        name=name,
        length=len(raw),
        start=(
            aligned_offset % len(template)
            if circular else aligned_offset
        ),
        end=(
            (aligned_offset + len(oriented)) % len(template) or len(template)
            if circular else min(len(template), aligned_offset + len(oriented))
        ),
        reverse_complemented=flipped,
        identity=identity,
        matched=matched,
        differences=differences,
        trimmed_start=cut_start,
        trimmed_end=cut_end,
        warnings=warnings,
        mean_quality=round(trace.mean_quality, 1) if trace is not None else None,
    )


def assemble_consensus(
    design: str,
    traces: dict[str, Chromatogram],
    *,
    circular: bool = False,
    region: tuple[int, int] | None = None,
    trim_quality: int = 0,
) -> ConsensusReport:

    template = clean_dna(design)
    if not template:
        raise SequenceError("The design is empty.")

    start, stop = region or (0, len(template))
    if not (0 <= start < stop <= len(template)):
        raise SequenceError("The consensus region is outside the design.")


    calls: dict[int, list[tuple[str, int, str, bool]]] = {
        position: [] for position in range(start, stop)
    }
    warnings: list[str] = []

    for name, trace in traces.items():
        raw = clean_dna(trace.sequence)
        if not raw:
            warnings.append(f"{name} is empty.")
            continue

        if trim_quality > 0 and trace.quality:
            cut_start, cut_stop = trace.trim(trim_quality)
        else:
            cut_start, cut_stop = 0, len(raw)
        trimmed = raw[cut_start:cut_stop]
        if len(trimmed) < ANCHOR:
            warnings.append(
                f"{name} has only {len(trimmed)} admitted bases and could not "
                "be placed in the consensus."
            )
            continue

        placed = _locate(template, trimmed, circular=circular)
        if placed is None:
            warnings.append(
                f"{name} does not match the reference closely enough to enter "
                "the consensus."
            )
            continue
        offset, flipped = placed
        full_oriented = reverse_complement(trimmed) if flipped else trimmed

        oriented_prefix = 0
        oriented_suffix = 0
        aligned_offset = offset
        oriented = full_oriented
        if not circular:
            oriented_prefix = max(0, -offset)
            oriented_suffix = max(0, offset + len(full_oriented) - len(template))
            oriented_stop = (
                len(full_oriented) - oriented_suffix
                if oriented_suffix else len(full_oriented)
            )
            oriented = full_oriented[oriented_prefix:oriented_stop]
            aligned_offset = max(0, offset)
            if len(oriented) < ANCHOR:
                warnings.append(
                    f"{name} overlaps only {len(oriented)} reference bases and "
                    "could not enter the consensus."
                )
                continue

        try:
            _operations, _matched, columns = _align(
                template, oriented, aligned_offset, circular=circular,
            )
        except SequenceError as error:
            warnings.append(f"{name}: {error}")
            continue
        for position, _expected, found, read_at in columns:
            if position is None or read_at is None or not found:
                continue
            if not start <= position < stop:
                continue
            full_oriented_at = oriented_prefix + read_at
            original_index = (
                len(full_oriented) - 1 - full_oriented_at
                if flipped else full_oriented_at
            ) + cut_start
            calls[position].append((
                found,
                trace.quality_at(original_index),
                name,
                flipped,
            ))

    positions: list[ConsensusPosition] = []
    consensus: list[str] = []
    missing: list[int] = []
    overlap_positions = 0
    agreeing_overlap_positions = 0
    identical_calls = 0
    called_positions = 0

    for position in range(start, stop):
        position_calls = calls[position]
        if not position_calls:
            consensus.append("N")
            missing.append(position)
            continue

        called_positions += 1
        by_base: dict[str, list[tuple[int, str]]] = {}
        for base, quality, name, _flipped in position_calls:
            by_base.setdefault(base, []).append((quality, name))


        ranked = sorted(
            (
                (len(support), sum(quality for quality, _name in support), base)
                for base, support in by_base.items()
            ),
            reverse=True,
        )
        best_count, best_quality, best_base = ranked[0]
        tied = (
            len(ranked) > 1
            and ranked[1][0] == best_count
            and ranked[1][1] == best_quality
        )
        call = "N" if tied else best_base
        agreement = len(by_base) == 1
        orientations = {flipped for _base, _quality, _name, flipped in position_calls}
        bidirectional = len(orientations) > 1
        if bidirectional:
            overlap_positions += 1
            if agreement:
                agreeing_overlap_positions += 1

        if call == template[position]:
            identical_calls += 1
        consensus.append(call)
        supporting = by_base.get(call, [])
        positions.append(ConsensusPosition(
            position=position,
            reference=template[position],
            call=call,
            supporting_reads=tuple(name for _quality, name in supporting),
            qualities=tuple(quality for quality, _name in supporting),
            combined_quality=min(60, best_quality) if call != "N" else 0,
            agreement=agreement,
        ))

    gaps: list[tuple[int, int]] = []
    for position in missing:
        if gaps and position == gaps[-1][1]:
            gaps[-1] = (gaps[-1][0], position + 1)
        else:
            gaps.append((position, position + 1))

    span = stop - start
    return ConsensusReport(
        reference_length=len(template),
        sequence="".join(consensus),
        positions=positions,
        gaps=gaps,
        coverage=round(100.0 * called_positions / span, 1),
        identity=round(100.0 * identical_calls / called_positions, 2)
        if called_positions else 0.0,
        bidirectional_overlap=round(100.0 * overlap_positions / span, 1),
        bidirectional_agreement=round(
            100.0 * agreeing_overlap_positions / overlap_positions, 2
        ) if overlap_positions else 0.0,
        warnings=warnings,
    )


def verify(
    design: str,
    reads: list[str] | dict[str, str],
    *,
    circular: bool = False,
    trim: int = 30,
    coding_start: int | None = None,
    coding_end: int | None = None,
    region: tuple[int, int] | None = None,
    traces: dict[str, Chromatogram] | None = None,
    trim_quality: int = DEFAULT_TRIM_QUALITY,
) -> VerificationReport:

    template = clean_dna(design)
    entries = reads if isinstance(reads, dict) else {
        f"read {i + 1}": sequence for i, sequence in enumerate(reads)
    }

    aligned: list[ReadAlignment] = []
    warnings: list[str] = []
    for name, sequence in entries.items():
        try:
            aligned.append(verify_read(
                template, sequence, name=name, circular=circular, trim=trim,
                coding_start=coding_start, coding_end=coding_end,
                trace=(traces or {}).get(name),
                trim_quality=trim_quality,
            ))
        except SequenceError as error:
            warnings.append(str(error))

    start, stop = region or (0, len(template))
    wanted = set(range(start, stop))
    for read in aligned:
        if read.end > read.start:
            wanted -= set(range(read.start, read.end))
        else:
            wanted -= set(range(read.start, len(template)))
            wanted -= set(range(read.end))

    gaps: list[tuple[int, int]] = []
    for position in sorted(wanted):
        if gaps and position == gaps[-1][1]:
            gaps[-1] = (gaps[-1][0], position + 1)
        else:
            gaps.append((position, position + 1))

    span = max(1, stop - start)
    return VerificationReport(
        design_length=len(template),
        reads=aligned,
        gaps=gaps,
        coverage=round(100.0 * (span - len(wanted)) / span, 1),
        warnings=warnings,
    )
