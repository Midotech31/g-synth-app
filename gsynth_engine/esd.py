from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

from gsynth_engine.constants import left_remainders
from gsynth_engine.constants import overhang as enzyme_overhang
from gsynth_engine.sequence import (
    SequenceError,
    gc_content,
    is_palindrome,
    longest_homopolymer,
    reverse_complement,
)
from gsynth_engine.ssd import SSDResult, design_small_sequence
from gsynth_engine.thermo import ANNEALING, melting_temperature

MIN_OVERHANG = 4
MAX_OVERHANG = 8


OVERHANG_SUPPLY: Final[dict[int, int]] = {4: 22, 5: 92, 6: 482}


def _describe_end(sequence: str, kind: str) -> str:

    return f"{sequence} ({kind} overhang)" if sequence else "blunt"


@dataclass(frozen=True)
class OligoPair:


    index: int
    name: str
    forward: str
    reverse: str
    top_start: int
    top_end: int
    left_overhang: str
    right_overhang: str
    is_first: bool
    is_last: bool
    bottom_offset: int = 0

    @property
    def left_overhang_strand(self) -> str:

        if self.bottom_offset > 0:
            return "top"
        return "bottom" if self.bottom_offset < 0 else "blunt"

    @property
    def right_overhang_strand(self) -> str:

        end_offset = (self.bottom_offset + len(self.reverse)) - len(self.forward)
        if end_offset > 0:
            return "bottom"
        return "top" if end_offset < 0 else "blunt"

    @property
    def forward_length(self) -> int:
        return len(self.forward)

    @property
    def reverse_length(self) -> int:
        return len(self.reverse)

    @property
    def forward_tm(self) -> float:

        return round(melting_temperature(self.forward, conditions=ANNEALING), 1)

    @property
    def reverse_tm(self) -> float:
        return round(melting_temperature(self.reverse, conditions=ANNEALING), 1)

    @property
    def duplex_gc(self) -> float:
        return round(gc_content(self.forward), 1)


@dataclass
class ESDResult:


    construct_forward: str
    construct_reverse: str
    fragments: list[OligoPair]
    overhang_length: int
    ssd: SSDResult
    warnings: list[str] = field(default_factory=list)

    @property
    def fragment_count(self) -> int:
        return len(self.fragments)

    @property
    def construct_length(self) -> int:
        return len(self.construct_forward)

    @property
    def oligo_count(self) -> int:
        return 2 * len(self.fragments)

    @property
    def junction_overhangs(self) -> list[str]:

        return [f.right_overhang for f in self.fragments[:-1]]

    @property
    def longest_oligo(self) -> int:
        return max(
            (max(len(f.forward), len(f.reverse)) for f in self.fragments),
            default=0,
        )

    @property
    def terminal_ends(self) -> tuple[tuple[str, str], tuple[str, str]]:

        if not self.fragments:
            return ("", "blunt"), ("", "blunt")

        top = "".join(f.forward for f in self.fragments)
        bottom = "".join(reverse_complement(f.reverse) for f in self.fragments)
        offset = self.fragments[0].bottom_offset


        if offset > 0:
            left = (top[:offset], "5'")
        elif offset < 0:
            left = (bottom[:-offset], "3'")
        else:
            left = ("", "blunt")

        top_end, bottom_end = len(top), offset + len(bottom)
        if bottom_end > top_end:
            right = (bottom[top_end - offset:], "5'")
        elif top_end > bottom_end:
            right = (top[bottom_end:], "3'")
        else:
            right = ("", "blunt")
        return left, right

    def verify(self) -> list[str]:

        problems: list[str] = []

        rebuilt_top = "".join(f.forward for f in self.fragments)
        if rebuilt_top != self.construct_forward:
            problems.append(
                "The forward oligos do not reassemble into the designed "
                f"construct ({len(rebuilt_top)} nt rebuilt vs "
                f"{len(self.construct_forward)} nt designed)."
            )


        rebuilt_bottom = "".join(
            reverse_complement(f.reverse) for f in self.fragments
        )
        expected_bottom = reverse_complement(self.construct_reverse)
        if rebuilt_bottom != expected_bottom:
            problems.append("The reverse oligos do not reassemble the bottom strand.")


        for left, right in zip(self.fragments, self.fragments[1:], strict=False):
            if left.right_overhang != right.left_overhang:
                problems.append(
                    f"Fragments {left.index} and {right.index} do not share a "
                    f"junction: {left.right_overhang!r} vs {right.left_overhang!r}."
                )


        (left_seq, left_kind), (right_seq, right_kind) = self.terminal_ends
        for side, enzyme, seen, kind in (
            ("left", self.ssd.left_enzyme, left_seq, left_kind),
            ("right", self.ssd.right_enzyme, right_seq, right_kind),
        ):
            wanted = enzyme_overhang(enzyme)
            if (seen, kind) != wanted:
                problems.append(
                    f"The {side}-hand end of the assembled fragments is "
                    f"{_describe_end(seen, kind)}, but {enzyme} leaves "
                    f"{_describe_end(*wanted)} — it would not ligate into a "
                    f"vector cut with {enzyme}."
                )


        seen: dict[str, int] = {}
        for fragment in self.fragments[:-1]:
            overhang = fragment.right_overhang
            if overhang in seen:
                problems.append(
                    f"Junctions {seen[overhang]} and {fragment.index} share the "
                    f"overhang {overhang} — fragments could ligate out of order."
                )
            seen[overhang] = fragment.index
        return problems


def _confusable_with(overhang: str) -> set[str]:

    ball: set[str] = set()
    for word in (overhang, reverse_complement(overhang)):
        ball.add(word)
        for i, base in enumerate(word):
            ball.update(word[:i] + other + word[i + 1:]
                        for other in "ACGT" if other != base)
    return ball


class _OverhangPool:


    def __init__(self, forbidden: set[str]) -> None:
        self.taken: set[str] = set()
        self.forbidden = set(forbidden)

        self._near_junction: set[str] = set()

        self._near_terminal: set[str] = set()
        for end in forbidden:
            self._near_terminal |= _confusable_with(end)

    def take(self, overhang: str) -> None:
        self.taken.add(overhang)
        self._near_junction |= _confusable_with(overhang)

    def problem(self, overhang: str) -> str | None:

        if is_palindrome(overhang):
            return "palindromic — it would anneal to itself"
        if overhang in self.taken or reverse_complement(overhang) in self.taken:
            return "already used at another junction"


        if overhang in self._near_junction:
            return "within one base of another junction — they could cross-ligate"
        if overhang in self.forbidden or reverse_complement(overhang) in self.forbidden:
            return "matches a terminal restriction overhang"
        if overhang in self._near_terminal:
            return "within one base of a terminal overhang — it could ligate into the vector"
        if longest_homopolymer(overhang) >= len(overhang):
            return "a homopolymer run"
        gc = sum(1 for base in overhang if base in "GC")
        if gc == 0:
            return "no G or C — the junction would be too weak"
        if gc == len(overhang):
            return "all G/C — prone to mispairing"
        return None


def _choose_junctions(
    top: str,
    *,
    count: int,
    overhang_length: int,
    ds_start: int,
    ds_end: int,
    forbidden: set[str],
    search_window: int,
) -> list[int]:

    if count <= 0:
        return []

    span = ds_end - ds_start
    junctions: list[int] = []
    pool = _OverhangPool(forbidden)

    for i in range(1, count + 1):
        ideal = ds_start + round(span * i / (count + 1))
        chosen: int | None = None
        last_reason = "no candidate position was examined"


        offsets = [0]
        for step in range(1, search_window + 1):
            offsets.extend((step, -step))

        for offset in offsets:
            position = ideal + offset


            if position <= ds_start or position + overhang_length >= ds_end:
                continue
            if junctions and position - junctions[-1] < overhang_length + 1:
                continue
            overhang = top[position : position + overhang_length]
            if len(overhang) < overhang_length:
                continue
            reason = pool.problem(overhang)
            if reason is None:
                chosen = position
                pool.take(overhang)
                break
            last_reason = reason

        if chosen is None:
            raise SequenceError(
                f"Could not place junction {i} of {count} near position {ideal}: "
                f"every candidate overhang within ±{search_window} nt was "
                f"unusable (last reason: {last_reason}). Try a different "
                f"overhang length, or a different fragment size."
            )
        junctions.append(chosen)

    return junctions


def _supply(overhang_length: int) -> int:

    if overhang_length in OVERHANG_SUPPLY:
        return OVERHANG_SUPPLY[overhang_length]

    return max(OVERHANG_SUPPLY.values()) * 4 ** (overhang_length - 6)


def _place_junctions(
    top: str,
    *,
    count: int,
    overhang_length: int,
    ds_start: int,
    ds_end: int,
    forbidden: set[str],
    search_window: int,
) -> tuple[list[int], int, list[str]]:

    notes: list[str] = []


    length = overhang_length
    while length < MAX_OVERHANG and _supply(length) < count:
        length += 1
    if length > overhang_length:
        notes.append(
            f"Overhangs widened to {length} nt: {count} junctions need more "
            f"distinct overhangs than {overhang_length} nt can supply "
            f"({_supply(overhang_length)}). The oligos stay the same length."
        )


    spacing = max(1, (ds_end - ds_start) // (count + 1))
    widest = max(search_window, spacing // 2)

    attempted: list[str] = []
    last: SequenceError | None = None
    for candidate_length in range(length, MAX_OVERHANG + 1):
        window = search_window
        while True:
            try:
                junctions = _choose_junctions(
                    top,
                    count=count,
                    overhang_length=candidate_length,
                    ds_start=ds_start,
                    ds_end=ds_end,
                    forbidden=forbidden,
                    search_window=window,
                )
            except SequenceError as exc:
                attempted.append(f"{candidate_length} nt within ±{window} nt")
                last = exc
            else:
                if window > search_window:
                    notes.append(
                        f"Junctions were allowed to move up to ±{window} nt from "
                        f"their ideal position (rather than ±{search_window}) to "
                        f"find usable overhangs, so fragment lengths vary more "
                        f"than usual."
                    )
                if candidate_length > length:
                    notes.append(
                        f"Overhangs widened to {candidate_length} nt: the "
                        f"sequence did not offer enough usable "
                        f"{length} nt overhangs."
                    )
                return junctions, candidate_length, notes

            if window >= widest:
                break
            window = min(widest, window * 2)


    region = top[ds_start:ds_end]
    words = {region[i : i + MAX_OVERHANG] for i in range(len(region) - MAX_OVERHANG)}
    consumed = 2 * (1 + 3 * MAX_OVERHANG)
    if len(words) < consumed * count:
        raise SequenceError(
            f"This sequence does not contain enough distinct subsequences to "
            f"order {count} junctions: across {len(region)} bases it offers "
            f"only {len(words)} different {MAX_OVERHANG} nt words, and each "
            f"junction placed rules out about {consumed} of them. Assemble it "
            "in blocks — design each half as its own construct and join them "
            "with the terminal enzymes — or use longer oligos, so there are "
            "fewer junctions to place."
        ) from last

    raise SequenceError(
        "Could not find a set of overhangs that assembles this sequence in one "
        f"order. Tried {', '.join(attempted)}. Try longer oligos, so there are "
        "fewer junctions to place."
    ) from last


def design_extended_sequence(
    sequence: str,
    *,
    enzyme_pair: str = "NdeI / XhoI",
    is_coding: bool = False,
    remove_stop: bool = False,
    cleavage_site: str | None = "Thrombin",
    include_his_tag: bool = True,
    include_linkers: bool = True,
    target_oligo_length: int = 90,
    overhang_length: int = 4,
    search_window: int = 12,
) -> ESDResult:

    if not MIN_OVERHANG <= overhang_length <= MAX_OVERHANG:
        raise SequenceError(
            f"Overhang length must be between {MIN_OVERHANG} and {MAX_OVERHANG} "
            f"nt for Extended Sequence Design — got {overhang_length}."
        )
    if target_oligo_length < 3 * overhang_length:
        raise SequenceError(
            f"Target oligo length ({target_oligo_length} nt) is too short for "
            f"{overhang_length} nt overhangs — use at least "
            f"{3 * overhang_length} nt."
        )

    ssd = design_small_sequence(
        sequence,
        enzyme_pair=enzyme_pair,
        is_coding=is_coding,
        remove_stop=remove_stop,
        cleavage_site=cleavage_site,
        include_his_tag=include_his_tag,
        include_linkers=include_linkers,
    )
    warnings = list(ssd.warnings)

    top = ssd.forward
    bottom = reverse_complement(ssd.reverse)


    fwd_remainder, rev_remainder = left_remainders(ssd.left_enzyme)
    offset = len(fwd_remainder) - len(rev_remainder)

    ds_start = max(0, offset)
    ds_end = min(len(top), offset + len(bottom))


    forbidden = {ssd.left_overhang, ssd.right_overhang}
    forbidden = {o for o in forbidden if o}


    fragment_count = max(1, round(len(top) / target_oligo_length))
    if fragment_count > 1:

        max_fragments = max(1, (ds_end - ds_start) // (3 * overhang_length))
        if fragment_count > max_fragments:
            fragment_count = max_fragments
            warnings.append(
                f"Fragment count reduced to {fragment_count} so each fragment "
                f"stays long enough to carry {overhang_length} nt overhangs."
            )

    junctions, overhang_length, notes = _place_junctions(
        top,
        count=fragment_count - 1,
        overhang_length=overhang_length,
        ds_start=ds_start,
        ds_end=ds_end,
        forbidden=forbidden,
        search_window=search_window,
    )
    warnings.extend(notes)


    top_cuts = [0, *junctions, len(top)]
    bottom_cuts = [offset, *[t + overhang_length for t in junctions], offset + len(bottom)]

    fragments: list[OligoPair] = []
    for i in range(len(top_cuts) - 1):
        t0, t1 = top_cuts[i], top_cuts[i + 1]
        b0, b1 = bottom_cuts[i], bottom_cuts[i + 1]

        forward_oligo = top[t0:t1]
        bottom_piece = bottom[b0 - offset : b1 - offset]
        reverse_oligo = reverse_complement(bottom_piece)

        is_first = i == 0
        is_last = i == len(top_cuts) - 2

        left_overhang = ssd.left_overhang if is_first else top[t0 : t0 + overhang_length]
        if is_last:
            right_overhang = ssd.right_overhang
        else:
            right_overhang = top[t1 : t1 + overhang_length]

        fragments.append(
            OligoPair(
                index=i + 1,
                name=f"F{i + 1}",
                forward=forward_oligo,
                reverse=reverse_oligo,
                top_start=t0,
                top_end=t1,
                left_overhang=left_overhang,
                right_overhang=right_overhang,
                is_first=is_first,
                is_last=is_last,
                bottom_offset=b0 - t0,
            )
        )

    plan = ESDResult(
        construct_forward=ssd.forward,
        construct_reverse=ssd.reverse,
        fragments=fragments,
        overhang_length=overhang_length,
        ssd=ssd,
        warnings=warnings,
    )


    problems = plan.verify()
    if problems:
        raise SequenceError(
            "Internal error — the designed fragments do not reassemble into "
            "the construct. Please report this with the input sequence. "
            + " ".join(problems)
        )

    longest = plan.longest_oligo
    if longest > 200:
        warnings.append(
            f"The longest oligo is {longest} nt. Most suppliers synthesise up "
            f"to ~200 nt reliably — consider a smaller target oligo length."
        )
    return plan
