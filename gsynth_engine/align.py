from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from gsynth_engine.sequence import SequenceError, clean_dna, reverse_complement

DATA = Path(__file__).parent / "data"


MAX_CELLS = 4_000_000


_SIMILAR_GROUPS = (
    "AGST",
    "ILMVF",
    "KRH",
    "DENQ",
    "FYW",
    "ST",
    "NQ",
)


@lru_cache(maxsize=2)
def blosum62() -> dict[str, dict[str, int]]:

    record = json.loads((DATA / "blosum62.json").read_text())
    return record["scores"]


@dataclass(frozen=True)
class Scoring:


    match: int = 5
    mismatch: int = -4
    gap_open: int = 10
    gap_extend: int = 1
    matrix: dict[str, dict[str, int]] | None = None

    def score(self, a: str, b: str) -> int:
        if self.matrix is not None:
            row = self.matrix.get(a)
            if row is not None and b in row:
                return row[b]
            return self.mismatch
        return self.match if a == b else self.mismatch


PROTEIN_SCORING = Scoring(match=1, mismatch=-4, gap_open=11, gap_extend=1)


@dataclass
class Alignment:


    top: str
    marks: str
    bottom: str
    score: int
    mode: str
    identities: int
    similarities: int
    gaps: int

    start_a: int
    end_a: int
    start_b: int
    end_b: int
    reverse_complemented: bool = False
    is_protein: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.top)

    @property
    def identity(self) -> float:

        return round(100.0 * self.identities / self.length, 1) if self.length else 0.0

    @property
    def similarity(self) -> float:

        if not self.length:
            return 0.0
        return round(100.0 * (self.identities + self.similarities) / self.length, 1)

    def rows(self, width: int = 60) -> list[dict[str, object]]:

        out: list[dict[str, object]] = []
        seen_a, seen_b = self.start_a, self.start_b

        for offset in range(0, self.length, width):
            top = self.top[offset : offset + width]
            bottom = self.bottom[offset : offset + width]
            a_bases = sum(1 for c in top if c != "-")
            b_bases = sum(1 for c in bottom if c != "-")

            out.append({
                "top": top,
                "marks": self.marks[offset : offset + width],
                "bottom": bottom,
                "top_start": seen_a + 1 if a_bases else None,
                "top_end": seen_a + a_bases,
                "bottom_start": seen_b + 1 if b_bases else None,
                "bottom_end": seen_b + b_bases,
            })
            seen_a += a_bases
            seen_b += b_bases
        return out

    def to_text(self, width: int = 60) -> str:

        pad = len(str(max(self.end_a, self.end_b))) + 1
        lines: list[str] = []
        for row in self.rows(width):
            label = f"{row['top_start']:>{pad}}" if row["top_start"] else " " * pad
            other = f"{row['bottom_start']:>{pad}}" if row["bottom_start"] else " " * pad
            lines.append(f"{label} {row['top']}  {row['top_end']}")
            lines.append(f"{' ' * pad} {row['marks']}")
            lines.append(f"{other} {row['bottom']}  {row['bottom_end']}")
            lines.append("")
        return "\n".join(lines).rstrip() + "\n"


def _similar(a: str, b: str) -> bool:
    return any(a in group and b in group for group in _SIMILAR_GROUPS)


def _clean(sequence: str, *, is_protein: bool, field_name: str) -> str:
    text = "".join(sequence.split()).upper() if is_protein else clean_dna(sequence)
    text = "".join(c for c in text if c.isalpha())
    if not text:
        raise SequenceError(f"The {field_name} is empty.")
    return text


def _gotoh(
    a: str, b: str, scoring: Scoring, mode: str,
) -> tuple[str, str, int, int, int, int, int]:

    n, m = len(a), len(b)
    open_cost, extend = scoring.gap_open, scoring.gap_extend
    NEG = -(1 << 30)

    local = mode == "local"
    free_start = mode in ("local", "semi-global")

    row_m = [NEG] * (m + 1)
    row_x = [NEG] * (m + 1)
    row_y = [NEG] * (m + 1)

    if mode == "global":
        row_m[0] = 0
        for j in range(1, m + 1):
            row_y[j] = -(open_cost + extend * j)
    else:

        for j in range(m + 1):
            row_m[j] = 0 if free_start else (0 if j == 0 else NEG)

    pointers: list[bytearray] = []
    best = (0 if free_start else NEG, 0, 0, 0)

    for i in range(1, n + 1):
        new_m = [NEG] * (m + 1)
        new_x = [NEG] * (m + 1)
        new_y = [NEG] * (m + 1)
        trace = bytearray(m + 1)

        if mode in ("global", "semi-global"):
            new_x[0] = -(open_cost + extend * i)
        elif free_start:
            new_m[0] = 0

        residue = a[i - 1]
        score_of = scoring.score

        for j in range(1, m + 1):

            up_left_m = row_m[j - 1]
            up_left_x = row_x[j - 1]
            up_left_y = row_y[j - 1]
            if up_left_m >= up_left_x and up_left_m >= up_left_y:
                diagonal, came_from = up_left_m, 0
            elif up_left_x >= up_left_y:
                diagonal, came_from = up_left_x, 1
            else:
                diagonal, came_from = up_left_y, 2

            value = diagonal + score_of(residue, b[j - 1])
            flags = came_from
            if local and value <= 0:
                value, flags = 0, 4
            new_m[j] = value


            opened = row_m[j] - open_cost - extend
            extended = row_x[j] - extend
            if extended > opened:
                new_x[j], _ = extended, None
                flags |= 8
            else:
                new_x[j] = opened


            opened_y = new_m[j - 1] - open_cost - extend
            extended_y = new_y[j - 1] - extend
            if extended_y > opened_y:
                new_y[j] = extended_y
                flags |= 16
            else:
                new_y[j] = opened_y

            trace[j] = flags

            if free_start and new_m[j] > best[0]:
                best = (new_m[j], i, j, 0)

        pointers.append(trace)
        row_m, row_x, row_y = new_m, new_x, new_y


    if mode == "global":
        end_i, end_j = n, m
        layer = max(range(3), key=lambda k: (row_m, row_x, row_y)[k][m])
        score = (row_m, row_x, row_y)[layer][m]
    elif mode == "semi-global":
        end_i = n
        layer, end_j = max(((k, j) for k in range(3) for j in range(m + 1)),
                          key=lambda cell: (row_m, row_x, row_y)[cell[0]][cell[1]])
        score = (row_m, row_x, row_y)[layer][end_j]
    else:
        score, end_i, end_j, layer = best


    top_out: list[str] = []
    bottom_out: list[str] = []
    i, j = end_i, end_j

    while i > 0 and j > 0:
        flags = pointers[i - 1][j]
        if layer == 0:
            if local and flags & 4:
                break
            top_out.append(a[i - 1])
            bottom_out.append(b[j - 1])
            layer = flags & 3
            i, j = i - 1, j - 1
        elif layer == 1:
            top_out.append(a[i - 1])
            bottom_out.append("-")
            layer = 1 if flags & 8 else 0
            i -= 1
        else:
            top_out.append("-")
            bottom_out.append(b[j - 1])
            layer = 2 if flags & 16 else 0
            j -= 1

    if mode == "global":
        while i > 0:
            top_out.append(a[i - 1])
            bottom_out.append("-")
            i -= 1
        while j > 0:
            top_out.append("-")
            bottom_out.append(b[j - 1])
            j -= 1
    elif mode == "semi-global":

        while i > 0:
            top_out.append(a[i - 1])
            bottom_out.append("-")
            i -= 1

    top_out.reverse()
    bottom_out.reverse()
    return "".join(top_out), "".join(bottom_out), score, i, end_i, j, end_j


def align(
    a: str,
    b: str,
    *,
    mode: str = "global",
    is_protein: bool = False,
    scoring: Scoring | None = None,
    try_reverse: bool = True,
) -> Alignment:

    if mode not in ("global", "local", "semi-global"):
        raise SequenceError(
            f"Unknown alignment mode {mode!r}. Use global, local or semi-global."
        )

    first = _clean(a, is_protein=is_protein, field_name="first sequence")
    second = _clean(b, is_protein=is_protein, field_name="second sequence")

    if len(first) * len(second) > MAX_CELLS:
        raise SequenceError(
            f"{len(first):,} × {len(second):,} is too large to align directly. "
            f"For long, nearly identical sequences use the verification tool, "
            f"which anchors first and aligns only the region that matters."
        )

    if scoring is None:
        scoring = (
            Scoring(
                matrix=blosum62(),
                match=PROTEIN_SCORING.match,
                mismatch=PROTEIN_SCORING.mismatch,
                gap_open=PROTEIN_SCORING.gap_open,
                gap_extend=PROTEIN_SCORING.gap_extend,
            )
            if is_protein else Scoring()
        )

    candidates = [(second, False)]
    if try_reverse and not is_protein:
        candidates.append((reverse_complement(second), True))

    best_result = None
    for candidate, flipped in candidates:
        result = _gotoh(first, candidate, scoring, mode)
        if best_result is None or result[2] > best_result[0][2]:
            best_result = (result, candidate, flipped)

    (top, bottom, score, start_a, end_a, start_b, end_b), used, flipped = best_result

    marks: list[str] = []
    identities = similarities = gaps = 0
    for x, y in zip(top, bottom, strict=False):
        if x == "-" or y == "-":
            marks.append(" ")
            gaps += 1
        elif x == y:
            marks.append("|")
            identities += 1
        elif is_protein and _similar(x, y):
            marks.append(":")
            similarities += 1
        else:
            marks.append(".")

    warnings: list[str] = []
    if flipped:
        warnings.append(
            "The second sequence aligns to the reverse complement of the "
            "first. It is the same sequence, read the other way round."
        )


    if flipped:
        start_b, end_b = len(used) - end_b, len(used) - start_b

    return Alignment(
        top=top, marks="".join(marks), bottom=bottom, score=score, mode=mode,
        identities=identities, similarities=similarities, gaps=gaps,
        start_a=start_a, end_a=end_a, start_b=start_b, end_b=end_b,
        reverse_complemented=flipped, is_protein=is_protein, warnings=warnings,
    )
