"""Sequencing primers — reading back what you built.

Once a clone is picked, it gets sequenced, and that needs primers that sit
outside the insert and read into it. The rules are not subtle but they are
easy to get wrong in a way that wastes a week:

* **Place it back from the target.** The first thirty to fifty bases after
  a sequencing primer are unusable. A primer that starts exactly at the
  insert reads the insert's beginning as noise — the one part everybody
  most wants to check, because that is where the ATG and the tag are.
* **It has to be unique.** A primer that binds twice in the plasmid gives a
  superimposed trace and no usable sequence. On a circular template that
  means checking both strands over the whole molecule, not just the region
  of interest.
* **Tm decides the run, not the guesswork.** Melting temperatures come from
  the nearest-neighbour model under sequencing conditions, the same engine
  the oligo design uses, so two numbers in the same project mean the same
  thing.

Long inserts need more than the two flanking primers: a Sanger read gives
perhaps 700 usable bases, so anything longer gets internal primers spaced
to overlap.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.sequence import (
    SequenceError,
    clean_dna,
    gc_content,
    longest_homopolymer,
    reverse_complement,
)
from gsynth_engine.thermo import BufferConditions, melting_temperature

#: A sequencing reaction: dilute primer, standard salt. Quoted with every
#: Tm so the number means something.
SEQUENCING = BufferConditions(
    name="sequencing reaction", oligo_nM=200.0, na_mM=50.0
)

#: How much sequence a Sanger read gives before it degrades, in practice.
READ_LENGTH = 700

#: Bases after the primer that come back as noise. A primer must sit at
#: least this far back from anything you want to read.
DEAD_ZONE = 50


@dataclass(frozen=True)
class Primer:
    """One sequencing primer, with everything needed to order and use it."""

    name: str
    sequence: str
    start: int              #: 0-based position of its 5' end in the template
    direction: int          #: 1 reads left→right, -1 reads right→left
    tm: float
    gc: float
    #: Where its usable read begins and ends in the template.
    reads_from: int
    reads_to: int
    warnings: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def as_row(self) -> dict[str, object]:
        return {
            "Name": self.name,
            "Sequence (5'->3')": self.sequence,
            "Length (nt)": self.length,
            "Tm (°C)": self.tm,
            "GC (%)": self.gc,
            "Direction": "forward" if self.direction == 1 else "reverse",
            "Reads": f"{self.reads_from + 1}–{self.reads_to}",
        }


@dataclass
class PrimerSet:
    """Every primer needed to read one region, and what they leave uncovered."""

    primers: list[Primer] = field(default_factory=list)
    target_start: int = 0
    target_end: int = 0
    gaps: list[tuple[int, int]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def covers_target(self) -> bool:
        return not self.gaps

    @property
    def as_rows(self) -> list[dict[str, object]]:
        return [primer.as_row for primer in self.primers]


def _is_unique(template: str, candidate: str, *, circular: bool) -> bool:
    """True when the primer binds in exactly one place, on either strand.

    A primer that binds twice gives a superimposed trace, which reads as a
    sequence full of double peaks and is usually blamed on the sample.
    """
    scan = template + template[: len(candidate) - 1] if circular else template
    hits = scan.count(candidate) + scan.count(reverse_complement(candidate))
    return hits == 1


def _score(candidate: str, tm: float, *, tm_target: float) -> float:
    """How good a primer this is. Lower is better."""
    gc = gc_content(candidate)
    penalty = abs(tm - tm_target) * 2.0
    penalty += abs(gc - 50.0) * 0.4
    # A G or C at the 3' end holds the primer down while the polymerase
    # starts; three or more in a row makes it stick too well elsewhere.
    if candidate[-1] not in "GC":
        penalty += 6.0
    if candidate[-4:].count("G") + candidate[-4:].count("C") >= 4:
        penalty += 4.0
    penalty += 3.0 * max(0, longest_homopolymer(candidate) - 4)
    return penalty


def _pick(
    template: str,
    *,
    anchor: int,
    direction: int,
    circular: bool,
    tm_min: float,
    tm_max: float,
    tm_target: float,
    length_range: tuple[int, int],
    search: int,
) -> tuple[str, int] | None:
    """Best primer whose 3' end sits near `anchor`, or None.

    `anchor` is where the primer should end reading *from*; the search walks
    outwards so a primer that cannot be made at the ideal spot is found a
    little further away rather than not at all.
    """
    length = len(template)
    best: tuple[float, str, int] | None = None

    for shift in range(0, search):
        for sign in ((0,) if shift == 0 else (-1, 1)):
            at = anchor + sign * shift
            for size in range(length_range[0], length_range[1] + 1):
                if direction == 1:
                    start = at - size
                else:
                    start = at

                if not circular and (start < 0 or start + size > length):
                    continue

                window = "".join(
                    template[(start + i) % length] for i in range(size)
                )
                candidate = window if direction == 1 else reverse_complement(window)
                if set(candidate) - set("ACGT"):
                    continue

                tm = melting_temperature(candidate, conditions=SEQUENCING)
                if not tm_min <= tm <= tm_max:
                    continue
                if not _is_unique(template, window, circular=circular):
                    continue

                penalty = _score(candidate, tm, tm_target=tm_target)
                if best is None or penalty < best[0]:
                    best = (penalty, candidate, start % length)
        if best is not None and shift >= 8:
            break

    return (best[1], best[2]) if best else None


def design_sequencing_primers(
    template: str,
    *,
    target_start: int,
    target_end: int,
    circular: bool = True,
    name: str = "seq",
    tm_min: float = 50.0,
    tm_max: float = 65.0,
    tm_target: float = 57.0,
    length_range: tuple[int, int] = (18, 26),
    margin: int = 80,
    read_length: int = READ_LENGTH,
    search: int = 60,
) -> PrimerSet:
    """Primers that read across a region, both strands, with internal ones.

    Args:
        target_start / target_end: the stretch that must be read — usually
            the insert.
        margin: how far back from the target each flanking primer sits. Must
            exceed the dead zone after a sequencing primer, or the start of
            the insert comes back as noise.

    Returns:
        A :class:`PrimerSet`. `gaps` lists any part of the target no primer
        reaches, so a long insert that needs another primer says so rather
        than looking finished.

    Raises:
        SequenceError: for a target that does not lie in the template.
    """
    sequence = clean_dna(template)
    length = len(sequence)
    if not sequence:
        raise SequenceError("The template is empty.")
    if not 0 <= target_start < target_end <= length:
        raise SequenceError(
            f"The target {target_start + 1}–{target_end} does not lie inside a "
            f"{length} bp template."
        )
    if margin < DEAD_ZONE:
        raise SequenceError(
            f"A margin of {margin} bp puts the primer inside the {DEAD_ZONE} bp "
            f"of noise that follows it — the start of the target would be "
            f"unreadable. Use at least {DEAD_ZONE}."
        )

    primers: list[Primer] = []
    warnings: list[str] = []

    def add(anchor: int, direction: int, label: str) -> None:
        found = _pick(
            sequence, anchor=anchor, direction=direction, circular=circular,
            tm_min=tm_min, tm_max=tm_max, tm_target=tm_target,
            length_range=length_range, search=search,
        )
        if found is None:
            warnings.append(
                f"No usable {label} primer near position {anchor + 1}: nothing "
                f"in that stretch is unique with a Tm of {tm_min:.0f}–"
                f"{tm_max:.0f} °C."
            )
            return
        candidate, start = found
        size = len(candidate)
        # A read runs round the origin as readily as any other stretch, so
        # the range wraps rather than clamping. Clamping quietly reported a
        # primer just past the origin as reading almost nothing.
        if direction == 1:
            first, last = start + size + DEAD_ZONE, start + size + read_length
        else:
            first, last = start - read_length, start - DEAD_ZONE

        if circular:
            reads_from, reads_to = first % length, last % length
        else:
            reads_from, reads_to = max(0, min(length, first)), max(0, min(length, last))

        primers.append(Primer(
            name=f"{name}_{label}",
            sequence=candidate,
            start=start,
            direction=direction,
            tm=round(melting_temperature(candidate, conditions=SEQUENCING), 1),
            gc=round(gc_content(candidate), 1),
            reads_from=reads_from,
            reads_to=reads_to,
        ))

    # The two flanking primers, each reading into the target from outside.
    add((target_start - margin) % length, 1, "F")
    add((target_end + margin) % length, -1, "R")

    # Internal primers when the target is longer than one read can cover.
    span = target_end - target_start
    usable = max(1, read_length - DEAD_ZONE)
    if span > usable:
        # Overlap by a quarter of a read, so no junction lands at an end.
        step = int(usable * 0.75)
        for index, position in enumerate(
            range(target_start + step, target_end - DEAD_ZONE, step), start=1
        ):
            add(position % length, 1, f"F{index + 1}")

    covered: set[int] = set()
    for primer in primers:
        if primer.reads_to > primer.reads_from:
            covered |= set(range(primer.reads_from, primer.reads_to))
        else:
            covered |= set(range(primer.reads_from, length))
            covered |= set(range(0, primer.reads_to))

    missing = sorted(set(range(target_start, target_end)) - covered)
    gaps: list[tuple[int, int]] = []
    for position in missing:
        if gaps and position == gaps[-1][1]:
            gaps[-1] = (gaps[-1][0], position + 1)
        else:
            gaps.append((position, position + 1))

    if gaps:
        warnings.append(
            f"{sum(b - a for a, b in gaps)} bp of the target is not reached by "
            f"any primer. Add one, or sequence from a nearer site."
        )

    return PrimerSet(
        primers=primers,
        target_start=target_start,
        target_end=target_end,
        gaps=gaps,
        warnings=warnings,
    )
