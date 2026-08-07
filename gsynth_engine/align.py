"""Pairwise alignment — comparing two sequences that are not the same thing.

Distinct from `verify`, which places a sequencing read against the construct
it is supposed to be. That one assumes near-identity and exploits it. This
one makes no such assumption: two enterocin genes from different strains, a
gene and the variant a supplier returned, a designed protein against its
homologue.

**Affine gaps.** One long gap is one evolutionary event; ten separate gaps
are ten. Scoring every gap position the same makes the algorithm prefer
scattered single-base gaps, which is biologically backwards and produces
alignments that look wrong to anyone who reads them. Opening a gap is
therefore expensive and extending it is cheap (Gotoh 1982).

**Three modes, because the question differs.**

* *global* — the whole of both, end to end. For two variants of one gene.
* *local* — the best-matching stretch, ignoring the rest. For finding a
  domain, or the one region two distant sequences share.
* *semi-global* — the whole of the shorter, placed in the longer without
  penalty at the ends. For asking where a gene sits in a plasmid, which
  global alignment answers badly and local answers incompletely.

**Protein scoring is read, not typed.** BLOSUM62 ships as data generated
from Biopython's copy of the published matrix. A substitution matrix
transcribed by hand is 576 chances to be quietly wrong in a way that shifts
every alignment slightly.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from gsynth_engine.sequence import SequenceError, clean_dna, reverse_complement

DATA = Path(__file__).parent / "data"

#: How much dynamic programming to allow. Two 1,500 nt genes is 2.25 M cells
#: and about a second; beyond a few million the wait stops being reasonable
#: and the caller is better served by a different tool.
MAX_CELLS = 4_000_000

#: Residues that substitute for one another readily enough to mark as
#: similar rather than different. Grouped by the property that matters.
_SIMILAR_GROUPS = (
    "AGST",      # small
    "ILMVF",     # hydrophobic
    "KRH",       # basic
    "DENQ",      # acidic and amide
    "FYW",       # aromatic
    "ST",        # hydroxyl
    "NQ",
)


@lru_cache(maxsize=2)
def blosum62() -> dict[str, dict[str, int]]:
    """The published matrix, as data rather than as a literal in this file."""
    record = json.loads((DATA / "blosum62.json").read_text())
    return record["scores"]


@dataclass(frozen=True)
class Scoring:
    """What a match, a mismatch and a gap are worth.

    Defaults are the ones EMBOSS uses for DNA, which is what most people's
    intuition about an identity percentage was formed on.
    """

    match: int = 5
    mismatch: int = -4
    gap_open: int = 10        #: cost of starting a gap, positive
    gap_extend: int = 1       #: cost of each further position, positive
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
    """Two sequences laid against each other, and how well they agree."""

    top: str                  #: sequence a, with gaps
    marks: str                #: '|' identical, ':' similar, ' ' neither
    bottom: str               #: sequence b, with gaps
    score: int
    mode: str
    identities: int
    similarities: int
    gaps: int
    #: Where the aligned stretch sits in each input, 0-based half-open.
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
        """Identities plus conservative substitutions. Protein only."""
        if not self.length:
            return 0.0
        return round(100.0 * (self.identities + self.similarities) / self.length, 1)

    def rows(self, width: int = 60) -> list[dict[str, object]]:
        """Wrapped for display, numbered in each sequence's own coordinates."""
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
        """The rendering people paste into a lab notebook."""
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
    """Affine-gap dynamic programming. Returns the two gapped strings, the
    score, and where the alignment starts and ends in each sequence.

    Three layers: `M` ends aligned, `X` ends in a gap in b, `Y` ends in a gap
    in a. The traceback has to remember **which layer it is in**, not merely
    which layer won at each cell. Reading the winner at every step leaves a
    gap the moment the aligned layer scores higher there, which fragments one
    twelve-base deletion into four — the exact thing affine penalties exist
    to prevent.

    So each cell keeps three decisions in one byte: where M came from, and
    whether X and Y were opened or extended. One byte per cell is what makes
    a few million cells affordable in Python.
    """
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
        # Free start along b: the alignment may begin anywhere in it.
        for j in range(m + 1):
            row_m[j] = 0 if free_start else (0 if j == 0 else NEG)

    pointers: list[bytearray] = []
    best = (0 if free_start else NEG, 0, 0, 0)      # score, i, j, layer

    for i in range(1, n + 1):
        new_m = [NEG] * (m + 1)
        new_x = [NEG] * (m + 1)
        new_y = [NEG] * (m + 1)
        trace = bytearray(m + 1)

        if mode == "global":
            new_x[0] = -(open_cost + extend * i)
        elif free_start:
            new_m[0] = 0

        residue = a[i - 1]
        score_of = scoring.score

        for j in range(1, m + 1):
            # ── M: the two residues are aligned ────────────────────────────
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
                value, flags = 0, 4                  # a fresh local start
            new_m[j] = value

            # ── X: a gap in b, running down ────────────────────────────────
            opened = row_m[j] - open_cost - extend
            extended = row_x[j] - extend
            if extended > opened:
                new_x[j], _ = extended, None
                flags |= 8
            else:
                new_x[j] = opened

            # ── Y: a gap in a, running across ──────────────────────────────
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

    # ── Where the alignment ends, and in which layer ───────────────────────
    if mode == "global":
        end_i, end_j = n, m
        layer = max(range(3), key=lambda k: (row_m, row_x, row_y)[k][m])
        score = (row_m, row_x, row_y)[layer][m]
    elif mode == "semi-global":
        end_i = n
        end_j = max(range(m + 1), key=lambda j: row_m[j])
        score, layer = row_m[end_j], 0
    else:
        score, end_i, end_j, layer = best

    # ── Traceback, staying inside a gap until it is genuinely over ─────────
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
        # The whole of `a` is used; whatever is left of it hangs off the end.
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
    """Align two sequences and report how well they agree.

    Args:
        mode: "global", "local" or "semi-global". See the module docstring —
            the right one depends on the question.
        try_reverse: for DNA, also try the reverse complement of `b` and keep
            whichever scores better. A gene cloned the other way round is not
            a different gene, and an alignment that misses that is useless.

    Raises:
        SequenceError: for empty input, an unknown mode, or a pair too large
            to align without an unreasonable wait.
    """
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

    # Report positions in the caller's own numbering, not the flipped copy.
    if flipped:
        start_b, end_b = len(used) - end_b, len(used) - start_b

    return Alignment(
        top=top, marks="".join(marks), bottom=bottom, score=score, mode=mode,
        identities=identities, similarities=similarities, gaps=gaps,
        start_a=start_a, end_a=end_a, start_b=start_b, end_b=end_b,
        reverse_complemented=flipped, is_protein=is_protein, warnings=warnings,
    )
