"""Primer pair design with Tm matching and cloning overhangs."""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.restriction.enzymes import get_enzyme
from gsynth_core.sequence.ops import ensure_dna, gc_content, reverse_complement
from gsynth_core.thermo.tm import melting_temperature


@dataclass(frozen=True, slots=True)
class PrimerParams:
    min_len: int = 18
    max_len: int = 30
    optimal_len: int = 22
    min_tm: float = 55.0
    max_tm: float = 65.0
    optimal_tm: float = 60.0
    max_tm_diff: float = 3.0
    min_gc: float = 40.0
    max_gc: float = 60.0
    max_homopolymer: int = 4
    require_gc_clamp_3prime: bool = True
    na_conc_mM: float = 50.0
    mg_conc_mM: float = 1.5
    dntp_conc_mM: float = 0.2
    dna_conc_nM: float = 500.0


@dataclass(frozen=True, slots=True)
class Primer:
    sequence: str
    tm: float
    gc: float
    length: int
    start: int
    end: int
    strand: str
    score: float


@dataclass(frozen=True, slots=True)
class PrimerPair:
    forward: Primer
    reverse: Primer
    product_size: int
    tm_difference: float
    score: float


def _homopolymer_ok(seq: str, max_run: int) -> bool:
    if len(seq) < 2:
        return True
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            run += 1
            if run > max_run:
                return False
        else:
            run = 1
    return True


def _gc_clamp_ok(seq: str) -> bool:
    n = sum(1 for c in seq[-5:] if c in "GC")
    return 1 <= n <= 3


def _score_primer(seq: str, tm: float, gc: float, p: PrimerParams) -> float:
    score = 100.0
    score -= 2.0 * (tm - p.optimal_tm) ** 2
    score -= 0.1 * (len(seq) - p.optimal_len) ** 2
    score -= 0.2 * abs(gc - 50.0)
    score += 5 if _gc_clamp_ok(seq) else -10
    if all(c in "AT" for c in seq[-3:]):
        score -= 15
    return score


def _design_one(template: str, *, strand: str, region: tuple[int, int],
                params: PrimerParams) -> Primer | None:
    best: Primer | None = None
    lo, hi = max(0, region[0]), min(len(template), region[1])
    for five in range(lo, hi):
        for L in range(params.min_len, params.max_len + 1):
            if strand == "+":
                if five + L > len(template):
                    continue
                seq = template[five:five+L]; ps, pe = five, five + L
            else:
                if five - L + 1 < 0:
                    continue
                seq = reverse_complement(template[five-L+1:five+1])
                ps, pe = five - L + 1, five + 1
            if not _homopolymer_ok(seq, params.max_homopolymer):
                continue
            g = gc_content(seq)
            if not (params.min_gc <= g <= params.max_gc):
                continue
            if params.require_gc_clamp_3prime and not _gc_clamp_ok(seq):
                continue
            tm = melting_temperature(seq, na_conc_mM=params.na_conc_mM,
                                     mg_conc_mM=params.mg_conc_mM,
                                     dntp_conc_mM=params.dntp_conc_mM,
                                     dna_conc_nM=params.dna_conc_nM)
            if not (params.min_tm <= tm <= params.max_tm):
                continue
            sc = _score_primer(seq, tm, g, params)
            if best is None or sc > best.score:
                best = Primer(sequence=seq, tm=tm, gc=g, length=L,
                              start=ps, end=pe, strand=strand, score=sc)
    return best


def design_primer_pair(template: str, *, region_to_amplify: tuple[int, int] | None = None,
                       params: PrimerParams | None = None,
                       anchor_window: int = 50) -> PrimerPair | None:
    """Design a Tm-matched forward/reverse primer pair."""
    template = ensure_dna(template)
    p = params or PrimerParams()
    L = len(template)
    lo_end, hi_end = (0, L) if region_to_amplify is None else region_to_amplify
    lo_end = max(0, lo_end); hi_end = min(L, hi_end)
    fwd = _design_one(template, strand="+",
                      region=(lo_end, min(L, lo_end + anchor_window)), params=p)
    rev = _design_one(template, strand="-",
                      region=(max(0, hi_end - anchor_window - 1), hi_end), params=p)
    if fwd is None or rev is None:
        return None
    tm_diff = abs(fwd.tm - rev.tm)
    if tm_diff > p.max_tm_diff:
        return None
    return PrimerPair(forward=fwd, reverse=rev,
                      product_size=rev.end - fwd.start,
                      tm_difference=tm_diff,
                      score=fwd.score + rev.score - 5 * tm_diff)


def design_cloning_primers(template: str, *, left_enzyme: str, right_enzyme: str,
                           region_to_amplify: tuple[int, int] | None = None,
                           params: PrimerParams | None = None,
                           flanking_nt: str = "AAAA") -> PrimerPair | None:
    """Append restriction sites + flanking nt for cloning."""
    p = params or PrimerParams()
    core = design_primer_pair(template, region_to_amplify=region_to_amplify, params=p)
    if core is None:
        return None
    le = get_enzyme(left_enzyme); re = get_enzyme(right_enzyme)
    fwd_seq = flanking_nt + le.site + core.forward.sequence
    rev_seq = flanking_nt + reverse_complement(re.site) + core.reverse.sequence
    from gsynth_core.thermo.tm import melting_temperature as tm_fn
    fwd_tm = tm_fn(core.forward.sequence, na_conc_mM=p.na_conc_mM, mg_conc_mM=p.mg_conc_mM,
                   dntp_conc_mM=p.dntp_conc_mM, dna_conc_nM=p.dna_conc_nM)
    rev_tm = tm_fn(core.reverse.sequence, na_conc_mM=p.na_conc_mM, mg_conc_mM=p.mg_conc_mM,
                   dntp_conc_mM=p.dntp_conc_mM, dna_conc_nM=p.dna_conc_nM)
    fwd = Primer(sequence=fwd_seq, tm=fwd_tm, gc=gc_content(fwd_seq),
                 length=len(fwd_seq), start=core.forward.start, end=core.forward.end,
                 strand="+", score=core.forward.score)
    rev = Primer(sequence=rev_seq, tm=rev_tm, gc=gc_content(rev_seq),
                 length=len(rev_seq), start=core.reverse.start, end=core.reverse.end,
                 strand="-", score=core.reverse.score)
    return PrimerPair(forward=fwd, reverse=rev,
                      product_size=core.product_size + len(fwd_seq) - len(core.forward.sequence)
                                                 + len(rev_seq) - len(core.reverse.sequence),
                      tm_difference=abs(fwd_tm - rev_tm), score=core.score)
