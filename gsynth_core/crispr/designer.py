"""Guide RNA design across 6 Cas variants with off-target scanning."""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from typing import Literal

from gsynth_core.crispr.scoring import cfd_score, mit_specificity_score, on_target_score
from gsynth_core.sequence.ops import find_motif_both_strands, gc_content, reverse_complement


class CasType(str, Enum):
    SpCas9 = "SpCas9"
    SpCas9_NG = "SpCas9-NG"
    SaCas9 = "SaCas9"
    Cas12a = "Cas12a"
    Cas13 = "Cas13"

    @property
    def pam(self) -> str:
        return {
            CasType.SpCas9: "NGG", CasType.SpCas9_NG: "NG",
            CasType.SaCas9: "NNGRRT", CasType.Cas12a: "TTTV", CasType.Cas13: "NNN",
        }[self]

    @property
    def protospacer_length(self) -> int:
        return {
            CasType.SpCas9: 20, CasType.SpCas9_NG: 20,
            CasType.SaCas9: 21, CasType.Cas12a: 23, CasType.Cas13: 28,
        }[self]

    @property
    def pam_downstream(self) -> bool:
        return self in (CasType.SpCas9, CasType.SpCas9_NG, CasType.SaCas9)


@dataclass(frozen=True, slots=True)
class Guide:
    protospacer: str
    pam: str
    start: int
    end: int
    strand: Literal["+", "-"]
    gc_percent: float
    on_target: float
    cas: CasType
    context: str
    specificity: float | None = None
    top_off_target_cfd: float | None = None
    off_target_count: int | None = None

    def __repr__(self) -> str:
        return (f"Guide({self.protospacer}+{self.pam}, {self.strand}{self.start}, "
                f"GC={self.gc_percent:.0f}%, on={self.on_target:.2f})")


def find_pam_sites(sequence: str, cas: CasType = CasType.SpCas9,
                   *, both_strands: bool = True) -> list[tuple[int, str]]:
    hits = find_motif_both_strands(sequence.upper(), cas.pam)
    out: list[tuple[int, str]] = [(p, "+") for p in hits["+"]]
    if both_strands:
        out.extend([(p, "-") for p in hits["-"]])
    return sorted(out)


def _extract_context(seq: str, ps: int, pe: int, strand: str, pam_len: int) -> str | None:
    """Extract a 30-mer Doench window: 4 nt upstream + 20 nt protospacer + PAM + tail.

    Doench 2016 was trained with a 3-nt PAM at positions 24..26 of the 30-mer.
    For variants with a shorter PAM (e.g. SpCas9-NG, 2 nt) we still return a
    30-nt window — the trailing downstream slice expands to keep the
    protospacer aligned to positions 4..23. Returns None if the window
    extends past either end of the input sequence.
    """
    if pam_len < 1:
        return None
    trail = max(0, 6 - pam_len)
    if strand == "+":
        s, e = ps - 4, pe + pam_len + trail
        if s < 0 or e > len(seq):
            return None
        return seq[s:e]
    rc = reverse_complement(seq)
    L = len(rc)
    rc_ps, rc_pe = L - pe, L - ps
    s, e = rc_ps - 4, rc_pe + pam_len + trail
    if s < 0 or e > L:
        return None
    return rc[s:e]


def design_guides(sequence: str, *, cas: CasType = CasType.SpCas9,
                  both_strands: bool = True, min_gc: float = 30.0, max_gc: float = 70.0,
                  min_on_target: float = 0.3, genome_reference: str | None = None,
                  max_off_targets: int = 50, top_n: int | None = None) -> list[Guide]:
    """Design and score gRNAs."""
    seq = sequence.upper()
    pam_len = len(cas.pam)
    pl = cas.protospacer_length
    candidates: list[Guide] = []
    for pam_pos, strand in find_pam_sites(seq, cas, both_strands=both_strands):
        if cas.pam_downstream:
            if strand == "+":
                ps, pe = pam_pos - pl, pam_pos
                if ps < 0 or pe > len(seq):
                    continue
                protospacer = seq[ps:pe]
                pam_seq = seq[pam_pos:pam_pos + pam_len]
            else:
                pfs = pam_pos + pam_len
                pfe = pfs + pl
                if pfe > len(seq):
                    continue
                protospacer = reverse_complement(seq[pfs:pfe])
                pam_seq = reverse_complement(seq[pam_pos:pam_pos + pam_len])
                ps, pe = pfs, pfe
        else:
            if strand == "+":
                ps = pam_pos + pam_len; pe = ps + pl
                if pe > len(seq):
                    continue
                protospacer = seq[ps:pe]
                pam_seq = seq[pam_pos:pam_pos + pam_len]
            else:
                pfe = pam_pos; pfs = pfe - pl
                if pfs < 0:
                    continue
                protospacer = reverse_complement(seq[pfs:pfe])
                pam_seq = reverse_complement(seq[pam_pos:pam_pos + pam_len])
                ps, pe = pfs, pfe
        gc = gc_content(protospacer)
        if not (min_gc <= gc <= max_gc):
            continue
        if cas in (CasType.SpCas9, CasType.SpCas9_NG):
            ctx = _extract_context(seq, ps, pe, strand, pam_len=pam_len)
            if ctx is None or len(ctx) != 30:
                continue
            on_score = on_target_score(ctx)
        else:
            ctx = ""
            on_score = 0.5 + 0.3 * (1 if 40 <= gc <= 60 else 0)
            if "TTTT" in protospacer:
                on_score *= 0.5
        if on_score < min_on_target:
            continue
        candidates.append(Guide(
            protospacer=protospacer, pam=pam_seq, start=ps, end=pe,
            strand=strand, gc_percent=gc, on_target=on_score, cas=cas, context=ctx,
        ))
    if genome_reference:
        candidates = _annotate_off_targets(candidates, genome_reference, max_off_targets)
    candidates.sort(key=lambda g: (-g.on_target, -(g.specificity or 0.0)))
    if top_n is not None:
        candidates = candidates[:top_n]
    return candidates


def _annotate_off_targets(guides: list[Guide], genome: str, max_per_guide: int) -> list[Guide]:
    gup = genome.upper()
    grc = reverse_complement(gup)
    out: list[Guide] = []
    for g in guides:
        proto = g.protospacer; L = len(proto)
        seed_len = 8
        seed = proto[-seed_len:]
        hits: list[tuple[str, str, int]] = []
        for ref in (gup, grc):
            idx = 0
            while True:
                pos = ref.find(seed, idx)
                if pos < 0:
                    break
                ps = pos - (L - seed_len)
                if ps < 0 or ps + L + len(g.pam) > len(ref):
                    idx = pos + 1; continue
                target = ref[ps:ps + L]
                tpam = ref[ps + L:ps + L + len(g.pam)]
                mm = sum(1 for a, b in zip(proto, target) if a != b)
                if mm == 0 and tpam == g.pam:
                    idx = pos + 1; continue
                hits.append((target, tpam, mm))
                if len(hits) >= max_per_guide:
                    break
                idx = pos + 1
            if len(hits) >= max_per_guide:
                break
        top_cfd = 0.0
        mm_pos_by_hit: list[list[int]] = []
        for target_seq, target_pam, _ in hits:
            if len(target_seq) == L:
                cfd = cfd_score(proto, target_seq, pam=target_pam if len(target_pam) == 3 else None)
                top_cfd = max(top_cfd, cfd)
                mm_pos_by_hit.append([i+1 for i, (a, b) in enumerate(zip(proto, target_seq)) if a != b])
        if mm_pos_by_hit:
            per_hit = [mit_specificity_score(m) / 100.0 for m in mm_pos_by_hit]
            spec = max(0.0, min(100.0, 100.0 / (1.0 + sum(per_hit))))
        else:
            spec = 100.0
        out.append(Guide(
            protospacer=g.protospacer, pam=g.pam, start=g.start, end=g.end,
            strand=g.strand, gc_percent=g.gc_percent, on_target=g.on_target,
            cas=g.cas, context=g.context, specificity=spec,
            top_off_target_cfd=top_cfd, off_target_count=len(hits),
        ))
    return out
