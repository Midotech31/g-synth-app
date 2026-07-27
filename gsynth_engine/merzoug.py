"""Merzoug Assembly — long constructs from annealed oligo pairs.

The method, as practised at the bench:

    Each fragment is a short duplex made by annealing one forward and one
    reverse oligo. Adjacent duplexes carry complementary 4–8 nt 5' overhangs,
    so they ligate in one defined order. Fragments are joined pairwise with
    ligase until the full-length insert is obtained, which is then cloned
    with the terminal restriction pair.

    No PCR at any step. No restriction enzyme is used on internal junctions —
    the overhangs are synthesised into the oligos directly.

**How the fragments are derived.** The full-length SSD duplex is built first
(see ssd.py), then *cut in silico*: the top strand is cut at a junction, the
bottom strand k nucleotides further along. That stagger is what leaves a
5' overhang on each side of the cut, and it guarantees two things for free —
the terminal ends are exactly the SSD enzyme overhangs, and re-ligating every
fragment reproduces the SSD construct base for base.

    top     5'──────────────┐ ┌──────────────3'
    bottom  3'──────────┐   └─┘          ┌───5'
                        └──── k nt ──────┘
                        the 5' overhang each junction presents

**Why junction placement is not arbitrary.** Overhangs decide the assembly
order. Two junctions sharing an overhang — or an overhang that is its own
reverse complement — let fragments ligate in the wrong order or to
themselves, and the error is invisible until sequencing. Junctions are
therefore nudged away from their ideal position until every overhang is
unique, non-palindromic, and distinct from the vector's own sticky ends.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from gsynth_engine.constants import left_remainders
from gsynth_engine.sequence import (
    SequenceError,
    gc_content,
    is_palindrome,
    longest_homopolymer,
    reverse_complement,
)
from gsynth_engine.thermo import ANNEALING, melting_temperature
from gsynth_engine.ssd import SSDResult, design_small_sequence

#: Overhang lengths the method allows. 4 nt is the practical minimum for a
#: stable ligation junction; beyond 8 nt the oligos get long without making
#: the junction more specific.
MIN_OVERHANG = 4
MAX_OVERHANG = 8


@dataclass(frozen=True)
class OligoPair:
    """One fragment: the two oligos to order, and how they anneal."""

    index: int                 # 1-based, the order fragments ligate in
    name: str
    forward: str               # top-strand oligo, 5'→3' — order this
    reverse: str               # bottom-strand oligo, 5'→3' — order this
    top_start: int             # position in the full construct (0-based)
    top_end: int
    left_overhang: str         # 5' overhang on the left, top-strand sense
    right_overhang: str        # 5' overhang on the right, top-strand sense
    is_first: bool
    is_last: bool

    @property
    def forward_length(self) -> int:
        return len(self.forward)

    @property
    def reverse_length(self) -> int:
        return len(self.reverse)

    @property
    def forward_tm(self) -> float:
        """Tm under the annealing reaction, not a generic primer dilution."""
        return round(melting_temperature(self.forward, conditions=ANNEALING), 1)

    @property
    def reverse_tm(self) -> float:
        return round(melting_temperature(self.reverse, conditions=ANNEALING), 1)

    @property
    def duplex_gc(self) -> float:
        return round(gc_content(self.forward), 1)


@dataclass
class AssemblyPlan:
    """Everything needed to order oligos and run the assembly."""

    construct_forward: str        # the full-length SSD forward strand
    construct_reverse: str        # the full-length SSD reverse strand
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
        """The internal overhangs, in assembly order."""
        return [f.right_overhang for f in self.fragments[:-1]]

    @property
    def longest_oligo(self) -> int:
        return max(
            (max(len(f.forward), len(f.reverse)) for f in self.fragments),
            default=0,
        )

    def verify(self) -> list[str]:
        """Re-derive the construct from the fragments. Returns problems found.

        An empty list means: ligating these oligos in order reproduces the
        designed construct exactly, on both strands.
        """
        problems: list[str] = []

        rebuilt_top = "".join(f.forward for f in self.fragments)
        if rebuilt_top != self.construct_forward:
            problems.append(
                "The forward oligos do not reassemble into the designed "
                f"construct ({len(rebuilt_top)} nt rebuilt vs "
                f"{len(self.construct_forward)} nt designed)."
            )

        # Bottom strand: read each reverse oligo back into top-strand sense.
        # Those pieces run left-to-right along the construct just as the
        # fragments do, so they concatenate in fragment order — the reverse
        # oligos are individually flipped, the series is not.
        rebuilt_bottom = "".join(
            reverse_complement(f.reverse) for f in self.fragments
        )
        expected_bottom = reverse_complement(self.construct_reverse)
        if rebuilt_bottom != expected_bottom:
            problems.append("The reverse oligos do not reassemble the bottom strand.")

        # Junctions must actually be complementary.
        for left, right in zip(self.fragments, self.fragments[1:]):
            if left.right_overhang != right.left_overhang:
                problems.append(
                    f"Fragments {left.index} and {right.index} do not share a "
                    f"junction: {left.right_overhang!r} vs {right.left_overhang!r}."
                )

        # Overhangs must be unique, or fragments ligate in the wrong order.
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


def _cross_ligates(overhang: str, other: str) -> bool:
    """True when two overhangs are close enough for T4 ligase to confuse them.

    NEB's ligase-fidelity work shows that overhangs differing at a single
    position still join at a measurable rate. Either strand can present the
    end, so both orientations are compared.
    """
    for candidate in (overhang, reverse_complement(overhang)):
        if len(candidate) != len(other):
            continue
        if sum(1 for a, b in zip(candidate, other) if a != b) <= 1:
            return True
    return False


def _overhang_problem(
    overhang: str, used: set[str], forbidden: set[str],
) -> str | None:
    """Why this overhang is unusable, or None when it is fine."""
    if is_palindrome(overhang):
        return "palindromic — it would anneal to itself"
    if overhang in used or reverse_complement(overhang) in used:
        return "already used at another junction"
    # A mispaired junction is a silent failure: the ligation works, the gel
    # looks right, and the error only appears at sequencing.
    if any(_cross_ligates(overhang, taken) for taken in used):
        return "within one base of another junction — they could cross-ligate"
    if overhang in forbidden or reverse_complement(overhang) in forbidden:
        return "matches a terminal restriction overhang"
    if any(_cross_ligates(overhang, end) for end in forbidden):
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
    """Pick `count` top-strand cut positions with well-behaved overhangs.

    Junctions start at evenly spaced ideal positions and are nudged outward,
    alternating left and right, until the overhang passes every rule.
    """
    if count <= 0:
        return []

    span = ds_end - ds_start
    junctions: list[int] = []
    used: set[str] = set()

    for i in range(1, count + 1):
        ideal = ds_start + round(span * i / (count + 1))
        chosen: int | None = None
        last_reason = "no candidate position was examined"

        # offsets: 0, +1, -1, +2, -2, …
        offsets = [0]
        for step in range(1, search_window + 1):
            offsets.extend((step, -step))

        for offset in offsets:
            position = ideal + offset
            # Both cuts must fall inside the double-stranded region, and
            # fragments must not collapse to nothing.
            if position <= ds_start or position + overhang_length >= ds_end:
                continue
            if junctions and position - junctions[-1] < overhang_length + 1:
                continue
            overhang = top[position : position + overhang_length]
            if len(overhang) < overhang_length:
                continue
            reason = _overhang_problem(overhang, used, forbidden)
            if reason is None:
                chosen = position
                used.add(overhang)
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


def design_merzoug_assembly(
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
) -> AssemblyPlan:
    """Design a full Merzoug assembly for one insert.

    Args:
        sequence: the insert. Everything about the cassette and the terminal
            sticky ends is handled by the SSD stage.
        target_oligo_length: how long each synthesised oligo should be. The
            real lengths land near this; junction placement moves them a
            little to keep the overhangs clean.
        overhang_length: 4–8 nt, per the method.
        search_window: how far a junction may move from its ideal position
            to find a usable overhang.

    Returns:
        An :class:`AssemblyPlan`. Call ``verify()`` before ordering — it
        re-ligates the fragments in silico and reports any mismatch.

    Raises:
        SequenceError: for input the user can fix, with a message saying how.
    """
    if not MIN_OVERHANG <= overhang_length <= MAX_OVERHANG:
        raise SequenceError(
            f"Overhang length must be between {MIN_OVERHANG} and {MAX_OVERHANG} "
            f"nt for Merzoug assembly — got {overhang_length}."
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
    bottom = reverse_complement(ssd.reverse)     # bottom strand, top-strand sense

    # Where the bottom strand sits relative to the top. The SSD construction
    # makes this the difference between the two left-hand enzyme remainders.
    fwd_remainder, rev_remainder = left_remainders(ssd.left_enzyme)
    offset = len(fwd_remainder) - len(rev_remainder)

    ds_start = max(0, offset)                      # first fully paired position
    ds_end = min(len(top), offset + len(bottom))   # one past the last

    # Terminal overhangs must never be reused internally, or a fragment could
    # ligate straight into the vector.
    forbidden = {ssd.left_overhang, ssd.right_overhang}
    forbidden = {o for o in forbidden if o}

    # How many fragments? Each forward oligo is roughly one fragment long.
    fragment_count = max(1, round(len(top) / target_oligo_length))
    if fragment_count > 1:
        # Keep every fragment long enough to carry two overhangs plus a core.
        max_fragments = max(1, (ds_end - ds_start) // (3 * overhang_length))
        if fragment_count > max_fragments:
            fragment_count = max_fragments
            warnings.append(
                f"Fragment count reduced to {fragment_count} so each fragment "
                f"stays long enough to carry {overhang_length} nt overhangs."
            )

    junctions = _choose_junctions(
        top,
        count=fragment_count - 1,
        overhang_length=overhang_length,
        ds_start=ds_start,
        ds_end=ds_end,
        forbidden=forbidden,
        search_window=search_window,
    )

    # Cut positions: top at t, bottom k further along.
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
            )
        )

    plan = AssemblyPlan(
        construct_forward=ssd.forward,
        construct_reverse=ssd.reverse,
        fragments=fragments,
        overhang_length=overhang_length,
        ssd=ssd,
        warnings=warnings,
    )

    # Never hand back a plan that does not re-ligate to the design.
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
