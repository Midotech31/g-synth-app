"""The hybridisation view — both strands, aligned, with the overhangs showing.

This is the check that has to happen before oligos are ordered. Every other
report is a summary; this one shows the actual molecule. It answers, at a
glance, the questions that cost a fortnight when they go unanswered:

* Do the two oligos of a fragment actually pair over their whole length?
* Does each junction present the 5' overhang it is supposed to, on the
  strand it is supposed to be on?
* Are the two terminal ends the vector's sticky ends, and nothing else?
* Where does one fragment stop and the next begin, relative to the ATG,
  the tag and the cleavage site?

**Coordinates.** One frame spans the whole construct. Column 0 is the 5'-most
base of the top strand — the first base of the left sticky end. The bottom
strand starts further right by the width of that sticky end, and runs further
right at the other end by the width of the right one. Wherever a strand is
absent the column holds a space, which is precisely what a single-stranded
overhang looks like::

    5'-TATGGGTTCTTCT…                 top strand starts here
          ||||||||||
    3'-  ACCCAAGAAGA…                 bottom starts two columns in

**Strand sense.** The top line reads 5'→3' left to right. The bottom line
reads 3'→5' left to right, which is how a duplex is drawn everywhere and how
the bases line up with their partners. The oligo you actually order for the
bottom strand is that line reversed — `OligoPair.reverse` holds it that way.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.merzoug import AssemblyPlan, OligoPair
from gsynth_engine.sequence import complement, reverse_complement

#: A column with no strand present.
GAP = " "

#: Line width for the text rendering. Sixty bases is the convention in
#: sequence viewers and fits an 80-column terminal with room for numbering.
DEFAULT_WIDTH = 60


@dataclass(frozen=True)
class Span:
    """A labelled stretch of the frame, half-open [start, end)."""

    name: str
    start: int
    end: int
    kind: str = "segment"

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class DuplexRow:
    """One wrapped line of the rendering."""

    start: int          #: frame coordinate of column 0 of this row
    top: str
    ticks: str          #: '|' where the strands pair
    bottom: str
    top_start: int | None    #: 1-based position of the first top base, if any
    top_end: int | None
    bottom_start: int | None
    bottom_end: int | None


@dataclass
class DuplexView:
    """Both strands of a construct in one coordinate frame."""

    top: str                              #: 5'→3', GAP where absent
    bottom: str                           #: 3'→5', GAP where absent
    segments: list[Span] = field(default_factory=list)
    top_fragments: list[Span] = field(default_factory=list)
    bottom_fragments: list[Span] = field(default_factory=list)
    junctions: list[int] = field(default_factory=list)
    left_overhang: str = ""
    right_overhang: str = ""

    @property
    def width(self) -> int:
        return len(self.top)

    def paired(self) -> str:
        """'|' where both strands are present and complementary, else ' '."""
        marks = []
        for top_base, bottom_base in zip(self.top, self.bottom):
            if top_base == GAP or bottom_base == GAP:
                marks.append(GAP)
            elif complement(top_base) == bottom_base:
                marks.append("|")
            else:
                # A mismatch here is a bug in the design, not in the drawing.
                marks.append("x")
        return "".join(marks)

    def mismatches(self) -> list[int]:
        """Frame coordinates where two present bases fail to pair."""
        return [i for i, mark in enumerate(self.paired()) if mark == "x"]

    def rows(self, width: int = DEFAULT_WIDTH) -> list[DuplexRow]:
        """Wrap the frame into fixed-width lines, numbered per strand.

        Numbering counts real bases, not columns, so a row that begins in the
        middle of an overhang still reports the position of its first base.
        """
        ticks = self.paired()
        rows: list[DuplexRow] = []

        for start in range(0, self.width, width):
            stop = min(start + width, self.width)
            top_chunk = self.top[start:stop]
            bottom_chunk = self.bottom[start:stop]

            rows.append(
                DuplexRow(
                    start=start,
                    top=top_chunk,
                    ticks=ticks[start:stop],
                    bottom=bottom_chunk,
                    top_start=_first_base_number(self.top, start, stop),
                    top_end=_last_base_number(self.top, start, stop),
                    bottom_start=_first_base_number(self.bottom, start, stop),
                    bottom_end=_last_base_number(self.bottom, start, stop),
                )
            )
        return rows

    def to_text(self, width: int = DEFAULT_WIDTH) -> str:
        """Plain-text rendering, for the protocol and for copy-paste."""
        number_width = len(str(self.width)) + 1
        lines: list[str] = []

        for row in self.rows(width):
            top_label = f"{row.top_start:>{number_width}}" if row.top_start else " " * number_width
            bottom_label = (
                f"{row.bottom_start:>{number_width}}" if row.bottom_start else " " * number_width
            )
            lines.append(f"{top_label} 5' {row.top}".rstrip())
            lines.append(f"{' ' * number_width}    {row.ticks}".rstrip())
            lines.append(f"{bottom_label} 3' {row.bottom}".rstrip())
            lines.append("")

        return "\n".join(lines).rstrip() + "\n"


def _first_base_number(strand: str, start: int, stop: int) -> int | None:
    """1-based index, among real bases, of the first base in [start, stop)."""
    before = sum(1 for ch in strand[:start] if ch != GAP)
    for offset in range(start, stop):
        if strand[offset] != GAP:
            return before + 1
    return None


def _last_base_number(strand: str, start: int, stop: int) -> int | None:
    count = sum(1 for ch in strand[:stop] if ch != GAP)
    return count if any(ch != GAP for ch in strand[start:stop]) else None


def _lay_out(top: str, bottom_sense: str, bottom_offset: int) -> tuple[str, str, int]:
    """Place two strands in one frame; returns the padded lines and top's column.

    `bottom_offset` is how far right the bottom strand's left end sits
    relative to the top strand's. It is positive for a 5' overhang on the
    left (the top strand starts first) and negative for a 3' overhang, where
    the bottom strand is the one that starts first — ApaI, KpnI, PstI and
    SacI all cut that way, and a view that assumed otherwise would draw
    their constructs wrong while still looking plausible.
    """
    top_column = max(0, -bottom_offset)
    bottom_column = max(0, bottom_offset)

    frame_width = max(top_column + len(top), bottom_column + len(bottom_sense))
    top_line = (GAP * top_column + top).ljust(frame_width, GAP)
    bottom_line = (GAP * bottom_column + complement(bottom_sense)).ljust(
        frame_width, GAP
    )
    return top_line, bottom_line, top_column


def construct_duplex(plan: AssemblyPlan) -> DuplexView:
    """Both strands of the full construct, with fragment boundaries marked.

    The two strands are staggered by the sticky ends the cloning enzymes
    leave, which is what makes the terminal overhangs and every internal
    junction visible as geometry rather than as a note in a table.
    """
    top = plan.construct_forward
    bottom_sense = reverse_complement(plan.construct_reverse)
    offset = plan.fragments[0].bottom_offset if plan.fragments else 0

    top_line, bottom_line, top_column = _lay_out(top, bottom_sense, offset)

    segments = [
        Span(
            name=segment.name,
            start=top_column + segment.start,
            end=top_column + segment.end,
        )
        for segment in plan.ssd.segments
    ]

    top_fragments = [
        Span(
            name=f.name,
            start=top_column + f.top_start,
            end=top_column + f.top_end,
            kind="fragment-top",
        )
        for f in plan.fragments
    ]

    # Walk the bottom pieces instead of deriving them from overhang widths:
    # they tile the bottom strand in fragment order, whatever the polarity of
    # the terminal ends.
    bottom_fragments: list[Span] = []
    cursor = max(0, offset)
    for fragment in plan.fragments:
        bottom_fragments.append(
            Span(
                name=fragment.name,
                start=cursor,
                end=cursor + len(fragment.reverse),
                kind="fragment-bottom",
            )
        )
        cursor += len(fragment.reverse)

    return DuplexView(
        top=top_line,
        bottom=bottom_line,
        segments=segments,
        top_fragments=top_fragments,
        bottom_fragments=bottom_fragments,
        junctions=[top_column + f.top_end for f in plan.fragments[:-1]],
        left_overhang=plan.ssd.left_overhang,
        right_overhang=plan.ssd.right_overhang,
    )


def fragment_duplex(fragment: OligoPair) -> DuplexView:
    """One fragment on its own, as it exists in its annealing tube.

    Built from the two oligos as ordered, not sliced out of the construct —
    so if the two ever disagreed, this view would show it.
    """
    top = fragment.forward
    bottom_sense = reverse_complement(fragment.reverse)
    top_line, bottom_line, top_column = _lay_out(
        top, bottom_sense, fragment.bottom_offset
    )

    return DuplexView(
        top=top_line,
        bottom=bottom_line,
        top_fragments=[
            Span(
                name=fragment.name,
                start=top_column,
                end=top_column + len(top),
                kind="fragment-top",
            )
        ],
        bottom_fragments=[
            Span(
                name=fragment.name,
                start=max(0, fragment.bottom_offset),
                end=max(0, fragment.bottom_offset) + len(bottom_sense),
                kind="fragment-bottom",
            )
        ],
        left_overhang=fragment.left_overhang,
        right_overhang=fragment.right_overhang,
    )
