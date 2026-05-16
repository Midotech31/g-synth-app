"""CRISPR sgRNA designer — pure-Python, no Streamlit/Reflex deps.

Same biology as `modules/crispr_designer.py` from the main G-Synth app,
condensed to what the POC needs:
- find SpCas9 NGG PAM sites on both strands
- extract 20-nt protospacers
- score by GC content + simple heuristics (TTTT terminator, GC clamp)
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Literal

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def clean_dna(seq: str) -> str:
    s = "".join(c for c in seq.upper() if c in "ACGT")
    return s


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in "GC")
    return 100.0 * gc / len(seq)


@dataclass
class Guide:
    protospacer: str
    pam: str
    start: int
    strand: Literal["+", "-"]
    gc_percent: float
    on_target: float
    has_tttt: bool

    def to_dict(self) -> dict:
        d = asdict(self)
        # round floats for compact JSON
        d["gc_percent"] = round(self.gc_percent, 1)
        d["on_target"] = round(self.on_target, 3)
        # Precompute UI badge color so Reflex foreach doesn't need to
        # do float-comparison on reactive vars (which is awkward in 0.9).
        if self.on_target >= 0.65:
            d["score_color"] = "green"
        elif self.on_target >= 0.50:
            d["score_color"] = "yellow"
        else:
            d["score_color"] = "red"
        d["strand_color"] = "green" if self.strand == "+" else "purple"
        return d


def _score(protospacer: str, gc: float) -> float:
    """Heuristic on-target score in [0, 1].

    Honest about being a heuristic — not a substitute for Doench 2016 /
    Azimuth in production, but enough to demonstrate the workflow.
    """
    score = 0.5
    # GC sweet spot 40–60%
    if 40.0 <= gc <= 60.0:
        score += 0.30
    elif 30.0 <= gc < 40.0 or 60.0 < gc <= 70.0:
        score += 0.10
    else:
        score -= 0.15
    # Penalize TTTT (Pol III terminator)
    if "TTTT" in protospacer:
        score -= 0.30
    # Light bonus if seed region (last 8 nt) has 4-6 GC
    seed = protospacer[-8:]
    seed_gc = sum(1 for c in seed if c in "GC")
    if 4 <= seed_gc <= 6:
        score += 0.05
    return max(0.0, min(1.0, score))


def design_guides(
    sequence: str,
    *,
    min_gc: float = 30.0,
    max_gc: float = 70.0,
    min_on_target: float = 0.4,
    top_n: int | None = 50,
) -> list[Guide]:
    """Find SpCas9 NGG-PAM guides on both strands, filter, sort, top-N."""
    seq = clean_dna(sequence)
    if len(seq) < 23:
        return []
    seqs = [(seq, "+"), (reverse_complement(seq), "-")]
    out: list[Guide] = []
    fwd_len = len(seq)
    for working, strand in seqs:
        L = len(working)
        for i in range(L - 22):
            window = working[i:i + 23]
            pam = window[20:23]
            # NGG match
            if pam[1] != "G" or pam[2] != "G":
                continue
            proto = window[:20]
            gc = gc_content(proto)
            if not (min_gc <= gc <= max_gc):
                continue
            score = _score(proto, gc)
            if score < min_on_target:
                continue
            # Coordinates in the forward strand
            start = i if strand == "+" else fwd_len - (i + 23)
            out.append(Guide(
                protospacer=proto,
                pam=pam,
                start=start,
                strand=strand,
                gc_percent=gc,
                on_target=score,
                has_tttt="TTTT" in proto,
            ))
    out.sort(key=lambda g: (-g.on_target, -g.gc_percent))
    if top_n is not None:
        out = out[:top_n]
    return out
