"""Sequencing verification — did you get the construct you designed?

The last step of the workflow, and the one still done by eye. A Sanger read
comes back, and someone squints at a chromatogram against a printout. This
compares the read to the designed sequence and says what differs, where,
and whether it matters to the protein.

**What makes it awkward.** A read arrives in whichever orientation the
primer happened to point, so half of them are the reverse complement of
what you are comparing against. The first thirty to fifty bases are noise
and the last two hundred degrade. And the construct is circular, so a read
across the origin is split in the coordinates but continuous in the
molecule. All three are handled here rather than left to the user, because
all three produce false differences that look exactly like real mutations.

**How it locates the read.** Exact k-mer anchors vote on an orientation and
an offset, then a banded alignment runs only in the located window. Full
dynamic programming over a 5 kb plasmid is 5 million cells for every read;
the band is a few thousand and gives the same answer, because a Sanger read
of a designed construct is nearly identical to it by construction.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from gsynth_engine.chromatogram import Chromatogram
from gsynth_engine.cloning import translate
from gsynth_engine.sequence import SequenceError, clean_dna, reverse_complement

#: Anchor length. Long enough to be unique in a plasmid, short enough that a
#: read with a mismatch every 40 bases still produces anchors.
ANCHOR = 14

#: How far the alignment may wander from the anchored diagonal. Sanger reads
#: of a designed construct do not carry large indels; anything bigger than
#: this is a different molecule, not a variant.
BAND = 24


@dataclass(frozen=True)
class Difference:
    """One place the read and the design disagree."""

    kind: str            #: "substitution", "insertion" or "deletion"
    position: int        #: 0-based, in the design
    expected: str
    found: str
    #: Filled in when the position falls inside a coding region.
    codon: int | None = None
    residue: int | None = None
    from_residue: str = ""
    to_residue: str = ""
    silent: bool | None = None
    #: Where this fell in the read, and how good that base was. Both are
    #: None when the read arrived as letters rather than as a trace.
    read_index: int | None = None
    quality: int | None = None

    @property
    def confident(self) -> bool | None:
        """Whether the trace supports this difference being real.

        Q20 is one wrong call in a hundred. Below it a "mutation" is more
        likely the basecaller choosing between a peak and its neighbour's
        shoulder, and the trace has to be looked at before anyone reorders
        an oligo over it.
        """
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
    """One sequencing read placed against the design."""

    name: str
    length: int
    #: Where the aligned part of the read sits in the design, 0-based.
    start: int
    end: int
    reverse_complemented: bool
    identity: float               #: percent over the aligned stretch
    matched: int
    differences: list[Difference] = field(default_factory=list)
    trimmed_start: int = 0        #: bases dropped from the read's 5' end
    trimmed_end: int = 0
    warnings: list[str] = field(default_factory=list)
    #: None when the read came as letters; a Phred mean when it came as a trace.
    mean_quality: float | None = None

    @property
    def unconfident_differences(self) -> list[Difference]:
        """Differences the trace does not support. Not mutations — noise."""
        return [d for d in self.differences if d.confident is False]

    @property
    def confirmed_differences(self) -> list[Difference]:
        """Differences a trace backs, plus every difference without one."""
        return [d for d in self.differences if d.confident is not False]

    @property
    def covered(self) -> int:
        return self.end - self.start

    @property
    def is_clean(self) -> bool:
        """True when the read agrees with the design everywhere it covers."""
        return not self.differences


@dataclass
class VerificationReport:
    """Every read, and what they say together."""

    design_length: int
    reads: list[ReadAlignment] = field(default_factory=list)
    #: Positions of the design no read covers.
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
        """Every read agrees, and something was actually read."""
        return bool(self.reads) and not self.differences

    @property
    def fully_covered(self) -> bool:
        return not self.gaps


def _anchors(sequence: str, k: int = ANCHOR) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i in range(len(sequence) - k + 1):
        index.setdefault(sequence[i : i + k], []).append(i)
    return index


def _locate(design: str, read: str, *, circular: bool) -> tuple[int, bool] | None:
    """Find the read's offset and orientation, or None when it does not fit.

    Anchors vote on a diagonal. A unique k-mer that lands at design position
    p from read position q implies an offset of p − q; the offset with the
    most votes wins. Repeated k-mers are ignored rather than allowed to vote
    for several diagonals at once.
    """
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


def _align(design: str, read: str, offset: int) -> tuple[list[tuple[int, str, str, int]], int]:
    """Banded global alignment of the read against the design at `offset`.

    Returns the edit operations as (design position, expected, found, read
    index) and how many positions matched. `expected` empty means an
    insertion in the read; `found` empty means a deletion. The read index is
    what lets a difference be shown against its own trace.

    **Everything here is sized by the band, not by the read.** Two rolling
    buffers, and a traceback row of the band's width indexed relative to its
    own start. Allocating a full row per base is what turns a linear
    algorithm back into a quadratic one: at 20,000 bases it was minutes of
    CPU and hundreds of megabytes for an answer the band computes in a
    fraction of a second.
    """
    n, m = len(read), len(read) + 2 * BAND
    window_start = offset - BAND
    length = len(design)
    window = "".join(design[(window_start + i) % length] for i in range(m))

    MATCH, MISMATCH, GAP = 1, -1, -2
    OUTSIDE = -(1 << 30)
    WIDTH = 2 * BAND + 4          # band, plus the cells its neighbours read

    previous = [j * GAP for j in range(m + 1)]
    current = [OUTSIDE] * (m + 1)
    trace: list[bytearray] = []
    bands: list[int] = []

    for i in range(1, n + 1):
        low = max(1, i - 1)
        high = min(m, i + 2 * BAND + 1)

        # Clear only what this row will write plus the cells the next row
        # reads past its end, so no value from an older row leaks in.
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
            best, step = diagonal, 0                   # 0 = match, 1 = ins, 2 = del
            if up > best:
                best, step = up, 1
            if left > best:
                best, step = left, 2
            current[j] = best
            row[j - low] = step

        trace.append(row)
        bands.append(low)
        previous, current = current, previous

    # Walk back from the best cell in the last row: the read is contained in
    # the window, so the design may extend past it on the right.
    last_low = bands[-1] if bands else 1
    last_high = min(m, n + 2 * BAND + 1)
    j = max(range(last_low, last_high + 1), key=lambda x: previous[x]) if n else 0
    i = n
    operations: list[tuple[int, str, str, int]] = []
    matched = 0

    while i > 0:
        low = bands[i - 1]
        if not low <= j < low + WIDTH:
            break                                      # fell out of the band
        step = trace[i - 1][j - low]
        if step == 0:
            expected, found = window[j - 1], read[i - 1]
            if expected == found:
                matched += 1
            else:
                operations.append(
                    ((window_start + j - 1) % length, expected, found, i - 1))
            i, j = i - 1, j - 1
        elif step == 1:
            operations.append(((window_start + j) % length, "", read[i - 1], i - 1))
            i -= 1
        else:
            operations.append(
                ((window_start + j - 1) % length, window[j - 1], "", max(0, i - 1)))
            j -= 1

    operations.reverse()
    return operations, matched


def _trim(read: str, *, trim: int) -> tuple[str, int, int]:
    """Drop the unreliable ends of a Sanger read.

    The first bases are the primer's own noise and the last few hundred
    degrade. Leaving them in produces differences that look exactly like
    real mutations, which is worse than reading less.
    """
    if trim <= 0 or len(read) <= 2 * trim:
        return read, 0, 0
    return read[trim : len(read) - trim], trim, trim


def _coding_effect(
    design: str, position: int, found: str, coding_start: int, coding_end: int,
) -> dict:
    """What a substitution does to the protein, when it lands in one."""
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
) -> ReadAlignment:
    """Place one sequencing read against the design and list what differs.

    Args:
        trim: bases to drop from each end before comparing. Sanger reads are
            noise at the start and degrade at the end; the default is
            deliberately conservative. Ignored when `trace` is given — the
            quality values say where the good sequence actually stops, and a
            fixed count is wrong in both directions on the same read.
        trace: the chromatogram this read came from. Supplying it trims by
            quality and marks each difference with the confidence of the base
            that produced it, which is what separates a mutation from a bad
            call.
        coding_start / coding_end: the reading frame, so a substitution can
            be reported as silent or as an amino-acid change.

    Raises:
        SequenceError: when the read cannot be placed in the design at all,
            which usually means it is a different construct.
    """
    template = clean_dna(design)
    raw = clean_dna(read)
    if not template:
        raise SequenceError("The design is empty.")
    if not raw:
        raise SequenceError(f"{name} is empty.")

    if trace is not None and trace.quality:
        # Mott trimming: the good stretch, not a fixed number of bases.
        start, stop = trace.trim()
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
    oriented = reverse_complement(trimmed) if flipped else trimmed

    operations, matched = _align(template, oriented, offset)

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
        index = len(oriented) - 1 - read_at if flipped else read_at
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
        start=offset % len(template),
        end=(offset + len(oriented)) % len(template) or len(template),
        reverse_complemented=flipped,
        identity=identity,
        matched=matched,
        differences=differences,
        trimmed_start=cut_start,
        trimmed_end=cut_end,
        warnings=warnings,
        mean_quality=round(trace.mean_quality, 1) if trace is not None else None,
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
) -> VerificationReport:
    """Check every read against the design and report coverage and changes.

    Args:
        region: the stretch that has to be covered — normally the insert.
            Coverage of a whole plasmid by two Sanger reads is never going
            to be complete, and reporting it as a failure would be noise.

    A read that cannot be placed becomes a warning rather than an exception:
    one bad trace out of six should not stop the other five being reported.
    """
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
            ))
        except SequenceError as error:
            warnings.append(str(error))

    start, stop = region or (0, len(template))
    wanted = set(range(start, stop))
    for read in aligned:
        if read.end > read.start:
            wanted -= set(range(read.start, read.end))
        else:                                    # wraps the origin
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
