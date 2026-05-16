"""Gibson Assembly with Tm-matched overlaps."""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.sequence.ops import ensure_dna
from gsynth_core.thermo.tm import melting_temperature


@dataclass(frozen=True, slots=True)
class GibsonFragment:
    name: str
    sequence: str
    left_overlap: str
    right_overlap: str
    left_tm: float
    right_tm: float
    order: int


def _find_overlap_length(left_seq: str, right_seq: str, *,
                        target_tm: float, tm_tolerance: float,
                        min_len: int, max_len: int) -> int:
    best_len = min_len
    best_diff = float("inf")
    for L in range(min_len, min(max_len, len(left_seq), len(right_seq)) + 1):
        tm = melting_temperature(left_seq[-L:])
        diff = abs(tm - target_tm)
        if diff < best_diff:
            best_diff = diff; best_len = L
        if tm > target_tm + tm_tolerance and L > min_len + 2:
            break
    return best_len


def design_gibson_overlaps(fragments: list[tuple[str, str]], *,
                          target_tm: float = 50.0, tm_tolerance: float = 5.0,
                          min_overlap_len: int = 15, max_overlap_len: int = 40,
                          circular: bool = True) -> list[GibsonFragment]:
    """Design Tm-matched overlaps for fragment series."""
    n = len(fragments)
    if n < 2:
        raise ValueError("Need at least 2 fragments")
    designed: list[GibsonFragment] = []
    for i, (name, seq) in enumerate(fragments):
        seq = ensure_dna(seq)
        prev_idx = i - 1 if i > 0 or circular else None
        next_idx = (i + 1) % n if i < n - 1 or circular else None
        if next_idx is not None:
            next_seq = ensure_dna(fragments[next_idx][1])
            L = _find_overlap_length(seq, next_seq, target_tm=target_tm,
                                    tm_tolerance=tm_tolerance,
                                    min_len=min_overlap_len, max_len=max_overlap_len)
            right_overlap = seq[-L:]
            right_tm = melting_temperature(right_overlap)
        else:
            right_overlap, right_tm = "", 0.0
        if prev_idx is not None:
            prev_seq = ensure_dna(fragments[prev_idx][1])
            L = _find_overlap_length(prev_seq, seq, target_tm=target_tm,
                                    tm_tolerance=tm_tolerance,
                                    min_len=min_overlap_len, max_len=max_overlap_len)
            left_overlap = prev_seq[-L:]
            left_tm = melting_temperature(left_overlap)
        else:
            left_overlap, left_tm = "", 0.0
        designed.append(GibsonFragment(
            name=name, sequence=seq,
            left_overlap=left_overlap, right_overlap=right_overlap,
            left_tm=left_tm, right_tm=right_tm, order=i,
        ))
    return designed


def simulate_gibson_assembly(fragments: list[GibsonFragment], *, circular: bool = True) -> str:
    """In-silico assembly: concatenate, folding shared overlaps."""
    if not fragments:
        return ""
    result = fragments[0].sequence
    for i in range(1, len(fragments)):
        f = fragments[i]
        if f.left_overlap and f.sequence.startswith(f.left_overlap):
            result += f.sequence[len(f.left_overlap):]
        else:
            result += f.sequence
    if circular and len(fragments) > 1:
        last_over = fragments[-1].right_overlap
        if last_over and result.endswith(last_over) and fragments[0].sequence.startswith(last_over):
            result = result[:-len(last_over)]
    return result
