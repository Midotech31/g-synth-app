"""Codon Adaptation Index (Sharp & Li 1987) and codon usage stats."""
from __future__ import annotations
import math
from collections import Counter
from typing import Mapping

from gsynth_core.constants import CODON_TABLES, GENETIC_CODE
from gsynth_core.sequence.ops import clean_dna


def relative_adaptiveness(organism: str = "e_coli") -> dict[str, float]:
    """w(codon) = f(codon) / f(max-codon for same aa) — for all 64 codons."""
    if organism not in CODON_TABLES:
        raise ValueError(f"Unknown organism {organism!r}")
    table = CODON_TABLES[organism]
    w: dict[str, float] = {}
    for aa, codons in table.items():
        fmax = max(codons.values())
        for codon, freq in codons.items():
            w[codon] = freq / fmax if fmax > 0 else 0.0
    return w


def calculate_cai(sequence: str, *, organism: str = "e_coli",
                  exclude_degenerate_one: bool = True) -> float:
    """CAI in (0, 1]. Excludes Met/Trp by default (Sharp & Li convention)."""
    w = relative_adaptiveness(organism)
    s = clean_dna(sequence)
    EPS = 0.01
    total, count = 0.0, 0
    for i in range(0, len(s) - 2, 3):
        codon = s[i:i+3]
        aa = GENETIC_CODE.get(codon)
        if aa is None or aa == "*":
            continue
        if exclude_degenerate_one and aa in ("M", "W"):
            continue
        total += math.log(max(w.get(codon, EPS), EPS))
        count += 1
    return math.exp(total / count) if count else 0.0


def codon_frequencies(sequence: str) -> dict[str, float]:
    """Empirical per-AA normalized codon frequencies."""
    s = clean_dna(sequence)
    counts: Counter[str] = Counter()
    for i in range(0, len(s) - 2, 3):
        counts[s[i:i+3]] += 1
    by_aa: dict[str, int] = {}
    for codon, n in counts.items():
        aa = GENETIC_CODE.get(codon)
        if aa:
            by_aa[aa] = by_aa.get(aa, 0) + n
    freqs: dict[str, float] = {}
    for codon, n in counts.items():
        aa = GENETIC_CODE.get(codon)
        if aa and by_aa.get(aa, 0) > 0:
            freqs[codon] = n / by_aa[aa]
    return freqs


def codon_usage_distance(sequence: str, *, organism: str = "e_coli") -> float:
    """RMS distance between observed codon usage and reference (0 = match, 1 = max)."""
    if organism not in CODON_TABLES:
        raise ValueError(f"Unknown organism {organism!r}")
    ref = CODON_TABLES[organism]
    obs = codon_frequencies(sequence)
    diffs: list[float] = []
    for aa, codons in ref.items():
        for codon, rf in codons.items():
            diffs.append((obs.get(codon, 0.0) - rf) ** 2)
    return math.sqrt(sum(diffs) / max(len(diffs), 1))


def codon_count_report(sequence: str, *, organism: str = "e_coli") -> list[dict]:
    """Per-codon report: count, observed/reference frequency, w."""
    if organism not in CODON_TABLES:
        raise ValueError(f"Unknown organism {organism!r}")
    ref = CODON_TABLES[organism]
    obs = codon_frequencies(sequence)
    w = relative_adaptiveness(organism)
    counts: Counter[str] = Counter()
    s = clean_dna(sequence)
    for i in range(0, len(s) - 2, 3):
        counts[s[i:i+3]] += 1
    report: list[dict] = []
    for aa, codons in ref.items():
        for codon, rf in codons.items():
            report.append({
                "codon": codon, "amino_acid": aa,
                "count": counts.get(codon, 0),
                "observed_freq": round(obs.get(codon, 0.0), 4),
                "reference_freq": rf,
                "w": round(w.get(codon, 0.0), 4),
            })
    report.sort(key=lambda r: (r["amino_acid"], -r["reference_freq"]))
    return report
