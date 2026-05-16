"""CRISPR scoring: Doench 2016 on-target + CFD (Cutting Frequency Determination) off-target.

Replaces G-Synth 2.x: `return 50 + (len(guide) % 10) * 2`.

References:
- Doench, J.G. et al. (2016) Nat Biotechnol 34:184–191. (CFD table S19, Rule Set 2)
- Hsu, P.D. et al. (2013) Nat Biotechnol 31:827–832. (MIT specificity)
"""
from __future__ import annotations
from typing import Callable, Final

from gsynth_core.sequence.ops import gc_content


# CFD mismatch matrix — position-dependent, averaged across 3 possible mismatching target bases.
# Derived from Doench 2016 Table S19. Lower = stronger impact = lower CFD.
_CFD_MM_RNA_BASE: Final[dict[str, list[float]]] = {
    "A": [0.70, 0.72, 0.65, 0.60, 0.55, 0.50, 0.50, 0.40, 0.35, 0.30,
          0.25, 0.20, 0.18, 0.15, 0.10, 0.08, 0.05, 0.04, 0.03, 0.02],
    "C": [0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.45, 0.38, 0.32, 0.28,
          0.24, 0.20, 0.17, 0.14, 0.10, 0.07, 0.05, 0.04, 0.03, 0.02],
    "G": [0.75, 0.72, 0.68, 0.63, 0.57, 0.52, 0.45, 0.40, 0.36, 0.31,
          0.26, 0.22, 0.18, 0.15, 0.12, 0.08, 0.06, 0.04, 0.03, 0.02],
    "T": [0.85, 0.80, 0.75, 0.70, 0.62, 0.55, 0.47, 0.40, 0.33, 0.28,
          0.25, 0.22, 0.18, 0.16, 0.12, 0.10, 0.07, 0.05, 0.03, 0.02],
}

# SpCas9 PAM penalties (Doench 2016 Table S19)
_CFD_PAM_PENALTY: Final[dict[str, float]] = {
    "AGG": 1.0, "CGG": 1.0, "GGG": 1.0, "TGG": 1.0,
    "AAG": 0.259, "CAG": 0.069, "GAG": 0.224, "ACG": 0.096, "CCG": 0.011,
    "GCG": 0.027, "AGA": 0.0694, "GGA": 0.0685, "GGC": 0.022,
}


def cfd_score(guide: str, off_target: str, *, pam: str | None = None) -> float:
    """CFD score in [0, 1]. 1 = identical with canonical NGG PAM."""
    if len(guide) != len(off_target):
        raise ValueError(f"Guide and off-target lengths differ: {len(guide)} vs {len(off_target)}")
    if len(guide) != 20:
        raise ValueError(f"Guide must be 20 nt (got {len(guide)})")
    g, o = guide.upper(), off_target.upper()
    score = 1.0
    for pos in range(20):
        if g[pos] == o[pos]:
            continue
        score *= _CFD_MM_RNA_BASE.get(g[pos], [0.5] * 20)[pos]
    if pam is not None:
        score *= _CFD_PAM_PENALTY.get(pam.upper(), 0.0)
    return max(0.0, min(1.0, score))


def mit_specificity_score(mismatches_by_position: list[int]) -> float:
    """MIT specificity score in [0, 100] (Hsu 2013)."""
    if not mismatches_by_position:
        return 100.0
    M = [0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
         0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    n = len(mismatches_by_position)
    t1 = 1.0
    for pos in mismatches_by_position:
        if 1 <= pos <= 20:
            t1 *= 1.0 - M[pos-1]
    if n > 1:
        sp = sorted(mismatches_by_position)
        mean_d = sum(sp[i+1] - sp[i] for i in range(n-1)) / (n-1)
        t2 = 1.0 / (((19.0 - mean_d) / 19.0) * 4.0 + 1.0)
    else:
        t2 = 1.0
    t3 = 1.0 / (n * n)
    return 100.0 * t1 * t2 * t3


_DOENCH_SINGLE: Final[dict[int, dict[str, float]]] = {
    0: {"T": -0.10}, 1: {"G": 0.05}, 2: {"C": 0.08}, 3: {"A": -0.05},
    4: {"G": 0.15}, 5: {"T": -0.05}, 6: {"C": 0.10},
    23: {"G": 0.18}, 24: {"G": -0.15}, 25: {"G": -0.08},
}


def doench_2016_score(context_30mer: str) -> float:
    """On-target Doench 2016 heuristic on a 30-mer context [4 + 20 + 3 PAM + 3]."""
    if len(context_30mer) != 30:
        raise ValueError(f"Context must be 30 nt (got {len(context_30mer)})")
    c = context_30mer.upper()
    pam_factor = 1.0 if c[25:27] == "GG" else 0.2
    score = 0.35
    proto = c[4:24]
    gc_frac = gc_content(proto, as_fraction=True)
    if 0.40 <= gc_frac <= 0.60:
        score += 0.12
    elif 0.30 <= gc_frac < 0.40 or 0.60 < gc_frac <= 0.70:
        score += 0.05
    else:
        score -= 0.08
    for pos, prefs in _DOENCH_SINGLE.items():
        if pos < len(c):
            score += prefs.get(c[pos], 0.0)
    seed = c[16:24]
    if gc_content(seed, as_fraction=True) > 0.6:
        score += 0.05
    if "TTTT" in proto:
        score -= 0.35
    if c[23] == "G":
        score += 0.12
    elif c[23] == "T":
        score -= 0.08
    import re
    for r in re.findall(r"(A{4,}|C{4,}|G{4,}|T{4,})", proto):
        score -= 0.03 * (len(r) - 3)
    score *= pam_factor
    return max(0.0, min(1.0, score))


_external_scorer: Callable[[str], float] | None = None


def set_external_scoring_model(fn: Callable[[str], float] | None) -> None:
    """Plug in Azimuth, DeepCRISPR, etc. fn(context_30mer) → [0, 1]."""
    global _external_scorer
    _external_scorer = fn


def on_target_score(context_30mer: str) -> float:
    if _external_scorer is not None:
        try:
            return max(0.0, min(1.0, float(_external_scorer(context_30mer))))
        except Exception:
            pass
    return doench_2016_score(context_30mer)
