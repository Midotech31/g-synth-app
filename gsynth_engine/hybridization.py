"""Antiparallel hybridization of two DNA strands.

Both inputs are written in the order supplied to an oligonucleotide vendor:
5′→3′.  The second strand is reverse-complemented before offsets are tested,
then drawn physically as 3′→5′ beneath the first.  This is deliberately not
pairwise homology alignment: internal gaps are not introduced, because a gap
would describe a bulge rather than a cohesive end.  Unpaired terminal bases
therefore remain visible as 5′ or 3′ overhangs.

The thermodynamic result is restricted to a perfectly complementary overlap.
SantaLucia nearest-neighbour parameters in :mod:`gsynth_engine.thermo` are not
silently applied to mismatched duplexes, for which this engine does not carry
the mismatch parameter set.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.sequence import (
    SequenceError,
    complement,
    reverse_complement,
    validate_dna,
)
from gsynth_engine.thermo import (
    ANNEALING,
    BufferConditions,
    duplex_thermodynamics,
    melting_temperature,
)

MAX_HYBRIDIZATION_CELLS = 4_000_000
MIN_COHESIVE_OVERLAP = 4


@dataclass(frozen=True)
class Overhang:
    """One single-stranded terminal extension, reported in strand sense."""

    end: str
    strand: str
    polarity: str
    sequence: str
    start: int
    end_position: int

    @property
    def length(self) -> int:
        return len(self.sequence)

    def payload(self) -> dict[str, object]:
        return {
            "end": self.end,
            "strand": self.strand,
            "polarity": self.polarity,
            "sequence": self.sequence,
            "length": self.length,
            "start": self.start,
            "end_position": self.end_position,
        }


@dataclass
class HybridizationResult:
    """A best ungapped antiparallel placement and its evidence."""

    first: str
    second: str
    top: str
    marks: str
    bottom: str
    offset: int
    overlap_start: int
    overlap_end: int
    paired_bases: int
    mismatches: int
    longest_perfect_run: int
    overhangs: list[Overhang]
    alternative_placements: int
    conditions: BufferConditions
    analysis_temperature_c: float
    tm_c: float | None = None
    tm_margin_c: float | None = None
    delta_h_kcal_mol: float | None = None
    delta_s_cal_mol_k: float | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def width(self) -> int:
        return len(self.top)

    @property
    def overlap_length(self) -> int:
        return max(0, self.overlap_end - self.overlap_start)

    @property
    def paired_percent(self) -> float:
        if not self.overlap_length:
            return 0.0
        return round(100.0 * self.paired_bases / self.overlap_length, 1)

    @property
    def complementarity(self) -> str:
        if self.mismatches == 0 and self.paired_bases >= MIN_COHESIVE_OVERLAP:
            return "exact"
        if self.paired_bases >= MIN_COHESIVE_OVERLAP:
            return "partial"
        return "insufficient"

    @property
    def predicted_state(self) -> str:
        if self.complementarity == "insufficient":
            return "insufficient_complementarity"
        if self.mismatches:
            return "mismatches_not_thermodynamically_scored"
        if self.tm_c is None:
            return "not_scored"
        if self.tm_c >= self.analysis_temperature_c:
            return "favourable_at_temperature"
        return "temperature_above_tm"

    @property
    def left_end(self) -> dict[str, object]:
        item = next(
            (overhang for overhang in self.overhangs if overhang.end == "left"),
            None,
        )
        return item.payload() if item else {"end": "left", "kind": "blunt", "length": 0}

    @property
    def right_end(self) -> dict[str, object]:
        item = next(
            (overhang for overhang in self.overhangs if overhang.end == "right"),
            None,
        )
        return item.payload() if item else {"end": "right", "kind": "blunt", "length": 0}

    def rows(self, width: int = 60) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        top_seen = bottom_seen = 0
        for start in range(0, self.width, width):
            stop = min(start + width, self.width)
            top = self.top[start:stop]
            bottom = self.bottom[start:stop]
            top_bases = sum(base != " " for base in top)
            bottom_bases = sum(base != " " for base in bottom)
            rows.append({
                "start": start,
                "stop": stop,
                "top": top,
                "marks": self.marks[start:stop],
                "bottom": bottom,
                "top_start": top_seen + 1 if top_bases else None,
                "top_end": top_seen + top_bases,
                # The lower strand is physically 3′→5′ in the drawing.
                "bottom_start": len(self.second) - bottom_seen if bottom_bases else None,
                "bottom_end": len(self.second) - bottom_seen - bottom_bases + 1,
            })
            top_seen += top_bases
            bottom_seen += bottom_bases
        return rows


def _longest_run(marks: str) -> int:
    best = current = 0
    for mark in marks:
        if mark == "|":
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def _candidate(
    first: str,
    second_rc: str,
    offset: int,
) -> tuple[tuple[int, ...], int, int, int]:
    """Return a sortable quality key and overlap measurements for one offset."""
    start = max(0, offset)
    stop = min(len(first), offset + len(second_rc))
    overlap = max(0, stop - start)
    paired = 0
    if overlap:
        paired = sum(
            first[index] == second_rc[index - offset]
            for index in range(start, stop)
        )
    mismatches = overlap - paired
    unpaired = len(first) + len(second_rc) - 2 * overlap
    score = 2 * paired - 3 * mismatches
    # Stable and deterministic: prefer more pairs, fewer mismatches, a longer
    # overlap, fewer dangling bases, then the placement closest to flush.
    key = (score, paired, -mismatches, overlap, -unpaired, -abs(offset), -offset)
    return key, start, stop, paired


def hybridize(
    first: str,
    second: str,
    *,
    conditions: BufferConditions = ANNEALING,
    analysis_temperature_c: float = 25.0,
) -> HybridizationResult:
    """Place two 5′→3′ DNA strands in their best antiparallel register.

    No internal gaps are introduced. Terminal displacement becomes explicit
    overhang geometry, while internal non-complementary columns remain marked
    as mismatches rather than being hidden by a local alignment.
    """
    top_sequence = validate_dna(first, field="first strand")
    partner_sequence = validate_dna(second, field="second strand")
    if len(top_sequence) * len(partner_sequence) > MAX_HYBRIDIZATION_CELLS:
        raise SequenceError(
            f"{len(top_sequence):,} × {len(partner_sequence):,} is too large for "
            "direct hybridization analysis. Analyse the intended oligonucleotide "
            "or cohesive-end region rather than an entire chromosome."
        )
    if not (-100.0 <= analysis_temperature_c <= 150.0):
        raise SequenceError("Analysis temperature must be between -100 and 150 °C.")

    second_rc = reverse_complement(partner_sequence)
    candidates: list[tuple[tuple[int, ...], int, int, int, int]] = []
    for offset in range(-len(second_rc) + 1, len(top_sequence)):
        key, start, stop, paired = _candidate(top_sequence, second_rc, offset)
        candidates.append((key, offset, start, stop, paired))
    candidates.sort(reverse=True)

    best_key, offset, overlap_start, overlap_end, paired = candidates[0]
    alternative_placements = (
        sum(candidate[0][:-2] == best_key[:-2] for candidate in candidates) - 1
    )

    union_start = min(0, offset)
    union_end = max(len(top_sequence), offset + len(second_rc))
    top_start = -union_start
    partner_start = offset - union_start
    width = union_end - union_start

    top = [" "] * width
    partner_rc = [" "] * width
    top[top_start : top_start + len(top_sequence)] = top_sequence
    partner_rc[partner_start : partner_start + len(second_rc)] = second_rc

    marks: list[str] = []
    bottom: list[str] = []
    for top_base, partner_base in zip(top, partner_rc, strict=True):
        bottom.append(complement(partner_base) if partner_base != " " else " ")
        if top_base == " " or partner_base == " ":
            marks.append(" ")
        elif top_base == partner_base:
            marks.append("|")
        else:
            marks.append("×")

    overlap_union_start = max(top_start, partner_start)
    overlap_union_end = min(
        top_start + len(top_sequence),
        partner_start + len(second_rc),
    )
    mismatch_count = sum(mark == "×" for mark in marks)

    overhangs: list[Overhang] = []
    if top_start < partner_start:
        sequence = "".join(top[top_start:partner_start])
        overhangs.append(
            Overhang("left", "first", "5′", sequence, top_start, partner_start)
        )
    elif partner_start < top_start:
        display = "".join(bottom[partner_start:top_start])
        overhangs.append(
            Overhang(
                "left",
                "second",
                "3′",
                display[::-1],
                partner_start,
                top_start,
            )
        )

    top_stop = top_start + len(top_sequence)
    partner_stop = partner_start + len(second_rc)
    if top_stop > partner_stop:
        sequence = "".join(top[partner_stop:top_stop])
        overhangs.append(
            Overhang("right", "first", "3′", sequence, partner_stop, top_stop)
        )
    elif partner_stop > top_stop:
        display = "".join(bottom[top_stop:partner_stop])
        overhangs.append(
            Overhang(
                "right",
                "second",
                "5′",
                display[::-1],
                top_stop,
                partner_stop,
            )
        )

    tm_c = margin = delta_h = delta_s = None
    warnings: list[str] = []
    if mismatch_count == 0 and paired:
        paired_sequence = "".join(top[overlap_union_start:overlap_union_end])
        tm_c = round(melting_temperature(paired_sequence, conditions=conditions), 1)
        margin = round(tm_c - analysis_temperature_c, 1)
        delta_h, delta_s = duplex_thermodynamics(paired_sequence)
        delta_h, delta_s = round(delta_h, 2), round(delta_s, 2)
    elif mismatch_count:
        warnings.append(
            "Tm is not reported for this mismatched overlap because the current "
            "nearest-neighbour model is parameterised for perfectly matched DNA."
        )

    if paired < MIN_COHESIVE_OVERLAP:
        warnings.append(
            f"Only {paired} complementary base pair(s) were found; fewer than "
            f"{MIN_COHESIVE_OVERLAP} bases is insufficient evidence for a useful cohesive end."
        )
    if alternative_placements:
        warnings.append(
            f"{alternative_placements + 1} placements have equivalent complementarity. "
            "Repeated sequence makes the register ambiguous."
        )
    if tm_c is not None and tm_c < analysis_temperature_c:
        warnings.append(
            f"The analysis temperature ({analysis_temperature_c:g} °C) is above "
            f"the predicted Tm ({tm_c:g} °C) under the stated conditions."
        )

    return HybridizationResult(
        first=top_sequence,
        second=partner_sequence,
        top="".join(top),
        marks="".join(marks),
        bottom="".join(bottom),
        offset=offset,
        overlap_start=overlap_union_start,
        overlap_end=overlap_union_end,
        paired_bases=paired,
        mismatches=mismatch_count,
        longest_perfect_run=_longest_run("".join(marks)),
        overhangs=overhangs,
        alternative_placements=alternative_placements,
        conditions=conditions,
        analysis_temperature_c=analysis_temperature_c,
        tm_c=tm_c,
        tm_margin_c=margin,
        delta_h_kcal_mol=delta_h,
        delta_s_cal_mol_k=delta_s,
        warnings=warnings,
    )
