"""Multi-objective codon optimization: CAI + GC window + homopolymer + repeat + restriction-site avoidance."""
from __future__ import annotations
import re
from dataclasses import dataclass, field

from gsynth_core.constants import CODON_TABLES, GENETIC_CODE, STANDARD_AA
from gsynth_core.codon.cai import calculate_cai
from gsynth_core.sequence.ops import gc_content, reverse_complement


@dataclass(frozen=True, slots=True)
class OptimizationParams:
    organism: str = "e_coli"
    gc_window_min: float = 40.0
    gc_window_max: float = 65.0
    gc_window_size: int = 30
    max_homopolymer: int = 5
    max_direct_repeat: int = 10
    avoid_sites: tuple[str, ...] = ()
    avoid_ribosome_binding: bool = True
    deterministic: bool = True
    weight_cai: float = 1.0
    weight_gc: float = 0.5
    weight_homopolymer: float = 2.0
    weight_repeat: float = 1.0


@dataclass(frozen=True, slots=True)
class OptimizationResult:
    protein: str
    original_dna: str | None
    optimized_dna: str
    cai_before: float | None
    cai_after: float
    gc_before: float | None
    gc_after: float
    changes: int
    forbidden_sites_remaining: list[str]
    warnings: list[str] = field(default_factory=list)


_IUPAC_MAP = {
    "A":"A","C":"C","G":"G","T":"T","R":"[AG]","Y":"[CT]","S":"[GC]","W":"[AT]",
    "K":"[GT]","M":"[AC]","B":"[CGT]","D":"[AGT]","H":"[ACT]","V":"[ACG]","N":"[ACGT]",
}


def _iupac_to_regex(pattern: str) -> str:
    return "".join(_IUPAC_MAP.get(c.upper(), re.escape(c)) for c in pattern)


def _site_in(seq: str, sites: tuple[str, ...]) -> list[str]:
    found: list[str] = []
    for site in sites:
        pat = _iupac_to_regex(site)
        rc = _iupac_to_regex(reverse_complement(site))
        if re.search(pat, seq) or re.search(rc, seq):
            found.append(site)
    return found


def _homopolymer_run(seq: str) -> int:
    if not seq:
        return 0
    last = seq[-1]
    run = 1
    for c in reversed(seq[:-1]):
        if c == last:
            run += 1
        else:
            break
    return run


def _has_direct_repeat(seq: str, min_len: int) -> bool:
    if len(seq) < 2 * min_len:
        return False
    for k in range(min_len, min(len(seq) // 2 + 1, 30)):
        if seq[-k:] == seq[-2*k:-k]:
            return True
    return False


def _score(codon: str, freq: float, context: str, p: OptimizationParams) -> tuple[float, list[str]]:
    candidate = context + codon
    violations: list[str] = []
    score = p.weight_cai * freq
    win = candidate[-p.gc_window_size:]
    if len(win) >= 15:
        gc = gc_content(win)
        if gc < p.gc_window_min:
            score -= p.weight_gc * ((p.gc_window_min - gc) / 100.0)
            if gc < p.gc_window_min - 15:
                violations.append(f"GC<{p.gc_window_min}")
        elif gc > p.gc_window_max:
            score -= p.weight_gc * ((gc - p.gc_window_max) / 100.0)
            if gc > p.gc_window_max + 15:
                violations.append(f"GC>{p.gc_window_max}")
    run = _homopolymer_run(candidate)
    if run > p.max_homopolymer:
        score -= p.weight_homopolymer * (run - p.max_homopolymer)
        violations.append(f"homopoly_{run}")
    if _has_direct_repeat(candidate, p.max_direct_repeat):
        score -= p.weight_repeat
        violations.append("repeat")
    if p.avoid_sites and _site_in(candidate[-20:], p.avoid_sites):
        score -= 10.0
        violations.append("restriction_site")
    if p.avoid_ribosome_binding and "AGGAGG" in candidate[-12:]:
        score -= 5.0
        violations.append("shine_dalgarno")
    return score, violations


def reverse_translate_best(protein: str, *, params: OptimizationParams | None = None) -> str:
    """Greedy best-in-context reverse translation (no post-pass)."""
    p = params or OptimizationParams()
    if p.organism not in CODON_TABLES:
        raise ValueError(f"Unknown organism {p.organism!r}")
    table = CODON_TABLES[p.organism]
    out_codons: list[str] = []
    context = ""
    for aa in protein.upper():
        if aa == "*":
            t = table["*"]; best = max(t, key=t.__getitem__)
            out_codons.append(best); context += best; continue
        if aa not in table:
            out_codons.append("NNN"); context += "NNN"; continue
        candidates = list(table[aa].items())
        if p.deterministic:
            candidates.sort(key=lambda x: (-x[1], x[0]))
        else:
            candidates.sort(key=lambda x: -x[1])
        best_score = float("-inf")
        best_codon = candidates[0][0]
        for codon, freq in candidates:
            s, _ = _score(codon, freq, context, p)
            if s > best_score:
                best_score = s; best_codon = codon
        out_codons.append(best_codon); context += best_codon
    return "".join(out_codons)


def _remove_forbidden(dna: str, protein: str, p: OptimizationParams,
                     max_passes: int = 5) -> tuple[str, list[str]]:
    if not p.avoid_sites:
        return dna, []
    table = CODON_TABLES[p.organism]
    current = dna
    for _ in range(max_passes):
        remaining = _site_in(current, p.avoid_sites)
        if not remaining:
            return current, []
        changed = False
        for site in list(set(remaining)):
            for pat in (site, reverse_complement(site)):
                m = re.search(_iupac_to_regex(pat), current)
                if not m:
                    continue
                lo, hi = m.span()
                first_c = lo // 3
                last_c = (hi - 1) // 3
                for cidx in range(first_c, last_c + 1):
                    cs = cidx * 3
                    if cs + 3 > len(current):
                        break
                    aa = protein[cidx] if cidx < len(protein) else None
                    if not aa or aa not in table:
                        continue
                    syns = [c for c in table[aa] if c != current[cs:cs+3]]
                    for alt in syns:
                        trial = current[:cs] + alt + current[cs+3:]
                        if not _site_in(trial[lo:hi+6], (site,)):
                            current = trial; changed = True; break
                    if changed:
                        break
                if changed:
                    break
            if changed:
                break
        if not changed:
            break
    return current, _site_in(current, p.avoid_sites)


def optimize(protein: str, *, params: OptimizationParams | None = None,
             original_dna: str | None = None) -> OptimizationResult:
    """Multi-objective codon optimization."""
    p = params or OptimizationParams()
    warnings: list[str] = []
    if not protein:
        raise ValueError("Empty protein sequence")
    invalid = set(protein.upper()) - set(STANDARD_AA + "*")
    if invalid:
        warnings.append(f"Non-standard AAs ignored/replaced with NNN: {sorted(invalid)}")
    opt_dna = reverse_translate_best(protein, params=p)
    opt_dna, remaining = _remove_forbidden(opt_dna, protein, p)
    if remaining:
        warnings.append(f"Could not remove all forbidden sites: {remaining}")
    cai_after = calculate_cai(opt_dna, organism=p.organism)
    gc_after = gc_content(opt_dna)
    cai_before = calculate_cai(original_dna, organism=p.organism) if original_dna else None
    gc_before = gc_content(original_dna) if original_dna else None
    changes = 0
    if original_dna and len(original_dna) == len(opt_dna):
        changes = sum(1 for a, b in zip(original_dna.upper(), opt_dna) if a != b)
    return OptimizationResult(
        protein=protein, original_dna=original_dna, optimized_dna=opt_dna,
        cai_before=cai_before, cai_after=cai_after,
        gc_before=gc_before, gc_after=gc_after,
        changes=changes, forbidden_sites_remaining=remaining, warnings=warnings,
    )
