"""Restriction digestion simulator."""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.restriction.enzymes import Enzyme, OverhangType, RestrictionSite, find_sites, get_enzyme


@dataclass(frozen=True, slots=True)
class DigestFragment:
    sequence: str
    length: int
    left_overhang: str
    right_overhang: str
    left_overhang_type: OverhangType
    right_overhang_type: OverhangType
    start_on_template: int
    end_on_template: int


def _cut_positions(sites: list[RestrictionSite]) -> list[tuple[int, int, Enzyme]]:
    seen: set[tuple[int, int]] = set()
    cuts: list[tuple[int, int, Enzyme]] = []
    for s in sites:
        key = (s.cut_position_fwd, s.cut_position_rev)
        if key not in seen:
            seen.add(key)
            cuts.append((s.cut_position_fwd, s.cut_position_rev, s.enzyme))
    cuts.sort(key=lambda x: (x[0], x[1]))
    return cuts


def _slice_overhang(seq: str, top: int, bot: int) -> tuple[str, OverhangType]:
    if top == bot:
        return "", OverhangType.BLUNT
    lo, hi = min(top, bot), max(top, bot)
    return seq[lo:hi], OverhangType.FIVE_PRIME if top < bot else OverhangType.THREE_PRIME


def _do_digest(seq: str, cuts: list[tuple[int, int, Enzyme]], circular: bool) -> list[DigestFragment]:
    L = len(seq)
    if not cuts:
        return [DigestFragment(
            sequence=seq, length=L, left_overhang="", right_overhang="",
            left_overhang_type=OverhangType.BLUNT, right_overhang_type=OverhangType.BLUNT,
            start_on_template=0, end_on_template=L,
        )]
    fragments: list[DigestFragment] = []
    if circular:
        first_top = cuts[0][0]
        rotated = seq[first_top:] + seq[:first_top]
        adj = sorted([((t - first_top) % L, (b - first_top) % L) for t, b, _ in cuts])
        ct = [c[0] for c in adj] + [L]
        cb = [c[1] for c in adj] + [L]
        for i in range(len(adj)):
            lt, lb = ct[i], cb[i]
            rt, rb = ct[i+1], cb[i+1]
            lo, lt_type = _slice_overhang(rotated, lt, lb)
            ro, rt_type = _slice_overhang(rotated, rt, rb)
            fragments.append(DigestFragment(
                sequence=rotated[lt:rt], length=rt - lt,
                left_overhang=lo, right_overhang=ro,
                left_overhang_type=lt_type, right_overhang_type=rt_type,
                start_on_template=(first_top + lt) % L,
                end_on_template=(first_top + rt) % L,
            ))
        return fragments
    boundaries = [(0, 0)] + [(t, b) for t, b, _ in cuts] + [(L, L)]
    for i in range(len(boundaries) - 1):
        at, ab = boundaries[i]
        bt, bb = boundaries[i+1]
        lo, lt_type = _slice_overhang(seq, at, ab)
        ro, rt_type = _slice_overhang(seq, bt, bb)
        fragments.append(DigestFragment(
            sequence=seq[at:bt], length=bt - at,
            left_overhang=lo, right_overhang=ro,
            left_overhang_type=lt_type, right_overhang_type=rt_type,
            start_on_template=at, end_on_template=bt,
        ))
    return fragments


def simulate_digestion(sequence: str, enzyme: str | Enzyme, *,
                       circular: bool = False) -> list[DigestFragment]:
    """Simulate single-enzyme digestion."""
    enz = enzyme if isinstance(enzyme, Enzyme) else get_enzyme(enzyme)
    sites = find_sites(sequence, [enz], both_strands=True)
    return _do_digest(sequence.upper(), _cut_positions(sites), circular)


def simulate_double_digestion(sequence: str, enzyme_a: str | Enzyme, enzyme_b: str | Enzyme,
                              *, circular: bool = False) -> list[DigestFragment]:
    """Simulate double digest (e.g. EcoRI + BamHI)."""
    ea = enzyme_a if isinstance(enzyme_a, Enzyme) else get_enzyme(enzyme_a)
    eb = enzyme_b if isinstance(enzyme_b, Enzyme) else get_enzyme(enzyme_b)
    sites = find_sites(sequence, [ea, eb], both_strands=True)
    return _do_digest(sequence.upper(), _cut_positions(sites), circular)
