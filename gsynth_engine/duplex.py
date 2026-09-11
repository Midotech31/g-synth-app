from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.esd import ESDResult, OligoPair
from gsynth_engine.sequence import SequenceError, complement, reverse_complement

GAP = " "


DEFAULT_WIDTH = 60


@dataclass(frozen=True)
class Span:


    name: str
    start: int
    end: int
    kind: str = "segment"

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class DuplexRow:


    start: int
    top: str
    ticks: str
    bottom: str
    top_start: int | None
    top_end: int | None
    bottom_start: int | None
    bottom_end: int | None


@dataclass
class DuplexView:


    top: str
    bottom: str
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

        marks = []
        for top_base, bottom_base in zip(self.top, self.bottom, strict=False):
            if top_base == GAP or bottom_base == GAP:
                marks.append(GAP)
            elif complement(top_base) == bottom_base:
                marks.append("|")
            else:

                marks.append("x")
        return "".join(marks)

    def mismatches(self) -> list[int]:

        return [i for i, mark in enumerate(self.paired()) if mark == "x"]

    def rows(self, width: int = DEFAULT_WIDTH) -> list[DuplexRow]:

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

    before = sum(1 for ch in strand[:start] if ch != GAP)
    for offset in range(start, stop):
        if strand[offset] != GAP:
            return before + 1
    return None


def _last_base_number(strand: str, start: int, stop: int) -> int | None:
    count = sum(1 for ch in strand[:stop] if ch != GAP)
    return count if any(ch != GAP for ch in strand[start:stop]) else None


def _lay_out(top: str, bottom_sense: str, bottom_offset: int) -> tuple[str, str, int]:

    top_column = max(0, -bottom_offset)
    bottom_column = max(0, bottom_offset)

    frame_width = max(top_column + len(top), bottom_column + len(bottom_sense))
    top_line = (GAP * top_column + top).ljust(frame_width, GAP)
    bottom_line = (GAP * bottom_column + complement(bottom_sense)).ljust(
        frame_width, GAP
    )
    return top_line, bottom_line, top_column


def construct_duplex(plan: ESDResult) -> DuplexView:

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


@dataclass
class JunctionView:


    name: str
    enzyme: str
    overhang: str
    kind: str
    compatible: bool
    reason: str = ""


    left_top: str = ""
    left_bottom: str = ""
    right_top: str = ""
    right_bottom: str = ""


    joined_top: str = ""
    joined_bottom: str = ""
    joined_pairs: str = ""
    seam: int = 0

    overhang_span: tuple[int, int] = (0, 0)

    @property
    def width(self) -> int:
        return len(self.joined_top)


def junction_view(
    plasmid: str,
    *,
    name: str,
    enzyme: str,
    position: int,
    overhang: str,
    kind: str,
    strand: str,
    flank: int = 18,
) -> JunctionView:

    length = len(plasmid)
    if length == 0:
        raise SequenceError("The plasmid is empty.")

    def at(index: int) -> str:
        return plasmid[index % length]

    width = len(overhang)
    joined_top = "".join(at(position + i) for i in range(-flank, flank))
    joined_bottom = complement(joined_top)
    seam = flank


    top_cut = seam
    if kind == "blunt":
        bottom_cut = seam
    elif kind == "5'":
        bottom_cut = seam + width
    else:
        bottom_cut = seam - width

    low, high = min(top_cut, bottom_cut), max(top_cut, bottom_cut)


    left_top = joined_top[:top_cut].ljust(high, GAP)
    left_bottom = joined_bottom[:bottom_cut].ljust(high, GAP)
    right_top = (GAP * (top_cut - low)) + joined_top[top_cut:]
    right_bottom = (GAP * (bottom_cut - low)) + joined_bottom[bottom_cut:]

    pairs = "".join(
        "|" if a != GAP and b != GAP and complement(a) == b else GAP
        for a, b in zip(joined_top, joined_bottom, strict=False)
    )

    span = (low, high)
    found = joined_top[low:high]
    compatible = kind == "blunt" or found == overhang
    reason = "" if compatible else (
        f"The plasmid reads {found} across the seam where {enzyme} leaves "
        f"{overhang}."
    )

    return JunctionView(
        name=name, enzyme=enzyme, overhang=overhang, kind=kind,
        compatible=compatible, reason=reason,
        left_top=left_top, left_bottom=left_bottom,
        right_top=right_top, right_bottom=right_bottom,
        joined_top=joined_top, joined_bottom=joined_bottom, joined_pairs=pairs,
        seam=seam, overhang_span=span,
    )
