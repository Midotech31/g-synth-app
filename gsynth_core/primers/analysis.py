"""Primer QC: hairpin, self-dimer, cross-dimer, GC clamp."""
from __future__ import annotations
from dataclasses import dataclass

from gsynth_core.sequence.ops import reverse_complement


@dataclass(frozen=True, slots=True)
class DimerRisk:
    max_complementary_run: int
    total_matches: int
    tm: float
    delta_g: float
    risk_level: str


def _try_primer3():
    try:
        import primer3
        return primer3
    except ImportError:
        return None


def _run_from_slides(a: str, b_revcomp: str) -> int:
    best = 0
    for shift in range(-len(b_revcomp) + 1, len(a)):
        run = max_run = 0
        for i in range(len(a)):
            j = i - shift
            if 0 <= j < len(b_revcomp):
                if a[i] == b_revcomp[j]:
                    run += 1; max_run = max(max_run, run)
                else:
                    run = 0
        best = max(best, max_run)
    return best


def _classify(max_run: int, tm: float, dg: float) -> str:
    if dg < -9.0 or tm > 45 or max_run >= 8:
        return "high"
    if dg < -5.0 or tm > 35 or max_run >= 5:
        return "medium"
    return "low"


def hairpin_risk(seq: str, *, tm_celsius: float = 60.0) -> DimerRisk:
    p3 = _try_primer3()
    s = seq.upper()
    if p3 is not None:
        try:
            hp = p3.calc_hairpin(s, temp_c=tm_celsius)
            tm = float(hp.tm) if hasattr(hp, "tm") and hp.tm else 0.0
            dg = float(hp.dg / 1000.0) if hasattr(hp, "dg") else 0.0
            return DimerRisk(0, 0, tm, dg, _classify(0, tm, dg))
        except Exception:
            pass
    mr = _run_from_slides(s, reverse_complement(s))
    return DimerRisk(mr, mr, 0.0, 0.0,
                     "high" if mr >= 6 else "medium" if mr >= 4 else "low")


def self_dimer_risk(seq: str, *, tm_celsius: float = 60.0) -> DimerRisk:
    p3 = _try_primer3()
    s = seq.upper()
    if p3 is not None:
        try:
            hd = p3.calc_homodimer(s, temp_c=tm_celsius)
            tm = float(hd.tm) if hd.tm else 0.0
            dg = float(hd.dg / 1000.0) if hd.dg else 0.0
            return DimerRisk(0, 0, tm, dg, _classify(0, tm, dg))
        except Exception:
            pass
    mr = _run_from_slides(s, reverse_complement(s))
    return DimerRisk(mr, mr, 0.0, 0.0,
                     "high" if mr >= 6 else "medium" if mr >= 4 else "low")


def cross_dimer_risk(fwd: str, rev: str, *, tm_celsius: float = 60.0) -> DimerRisk:
    p3 = _try_primer3()
    if p3 is not None:
        try:
            hd = p3.calc_heterodimer(fwd.upper(), rev.upper(), temp_c=tm_celsius)
            tm = float(hd.tm) if hd.tm else 0.0
            dg = float(hd.dg / 1000.0) if hd.dg else 0.0
            return DimerRisk(0, 0, tm, dg, _classify(0, tm, dg))
        except Exception:
            pass
    mr = _run_from_slides(fwd.upper(), reverse_complement(rev.upper()))
    return DimerRisk(mr, mr, 0.0, 0.0,
                     "high" if mr >= 6 else "medium" if mr >= 4 else "low")


@dataclass(frozen=True, slots=True)
class PrimerAnalysis:
    sequence: str
    length: int
    tm: float
    gc_content: float
    hairpin: DimerRisk
    self_dimer: DimerRisk
    gc_clamp_3prime: bool
    homopolymer_max: int


def analyze_primer(seq: str, *, annealing_temp: float = 60.0) -> PrimerAnalysis:
    """Full QC report on one primer."""
    from gsynth_core.thermo.tm import melting_temperature
    from gsynth_core.sequence.ops import gc_content
    s = seq.upper()
    last5 = s[-5:] if len(s) >= 5 else s
    gc3 = sum(1 for c in last5 if c in "GC")
    clamp_ok = 1 <= gc3 <= 3
    max_run = 1; run = 1
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            run += 1; max_run = max(max_run, run)
        else:
            run = 1
    return PrimerAnalysis(
        sequence=s, length=len(s),
        tm=melting_temperature(s), gc_content=gc_content(s),
        hairpin=hairpin_risk(s, tm_celsius=annealing_temp),
        self_dimer=self_dimer_risk(s, tm_celsius=annealing_temp),
        gc_clamp_3prime=clamp_ok, homopolymer_max=max_run,
    )
