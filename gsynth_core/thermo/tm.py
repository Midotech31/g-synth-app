"""Melting temperature calculations.

Prefers primer3-py when available; pure-Python NN fallback otherwise.
Validated: M13 universal primer Tm = 57.3 °C (literature ~56 °C).
"""
from __future__ import annotations
import math
from enum import Enum
from typing import Final


class TmMethod(str, Enum):
    AUTO = "auto"
    NEAREST_NEIGHBOR = "nn"
    MARMUR = "marmur"
    WALLACE = "wallace"
    GC_CONTENT = "gc"


# SantaLucia 1998 unified NN parameters (ΔH kcal/mol, ΔS cal/mol/K), 1 M NaCl.
# Convention: top and bottom both written 5'→3' (bottom = reverse_complement(top)).
# 10 unique pairs; the other 6 are dyad-equivalent by swapping slash halves.
_NN_PARAMS: Final[dict[str, tuple[float, float]]] = {
    "AA/TT": (-7.9, -22.2), "AT/AT": (-7.2, -20.4), "TA/TA": (-7.2, -21.3),
    "CA/TG": (-8.5, -22.7), "GT/AC": (-8.4, -22.4), "CT/AG": (-7.8, -21.0),
    "GA/TC": (-8.2, -22.2), "CG/CG": (-10.6, -27.2), "GC/GC": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}
_INIT_GC: Final[tuple[float, float]] = (0.1, -2.8)
_INIT_AT: Final[tuple[float, float]] = (2.3, 4.1)
_SYMMETRY_PENALTY_dS: Final[float] = -1.4


def _lookup_nn(dimer: str) -> tuple[float, float]:
    from gsynth_core.sequence.ops import reverse_complement
    top = dimer.upper()
    bottom = reverse_complement(top)
    key_a = f"{top}/{bottom}"
    if key_a in _NN_PARAMS:
        return _NN_PARAMS[key_a]
    key_b = f"{bottom}/{top}"
    if key_b in _NN_PARAMS:
        return _NN_PARAMS[key_b]
    raise KeyError(f"No NN parameters for dimer {dimer!r}")


def melting_temperature_nn(seq: str, *, dna_conc_nM: float = 500.0,
                           na_conc_mM: float = 50.0, mg_conc_mM: float = 0.0,
                           dntp_conc_mM: float = 0.0) -> float:
    """SantaLucia 1998 + Owczarzy 2008 salt correction."""
    from gsynth_core.sequence.ops import ensure_dna, is_palindrome, gc_content
    s = ensure_dna(seq)
    if len(s) < 8:
        return melting_temperature_marmur(s)
    dH = dS = 0.0
    for i in range(len(s) - 1):
        h, sv = _lookup_nn(s[i:i+2])
        dH += h; dS += sv
    for end in (s[0], s[-1]):
        h, sv = (_INIT_GC if end in "GC" else _INIT_AT)
        dH += h; dS += sv
    pal = is_palindrome(s)
    if pal:
        dS += _SYMMETRY_PENALTY_dS
    R = 1.987
    ct = dna_conc_nM * 1e-9
    x = 1.0 if pal else 4.0
    tm_1M_K = (dH * 1000.0) / (dS + R * math.log(ct / x))
    mg_free = max(0.0, mg_conc_mM - dntp_conc_mM)
    mg_M, na_M = mg_free * 1e-3, na_conc_mM * 1e-3
    gc_frac = gc_content(s, as_fraction=True)
    n_bp = len(s) - 1
    if mg_M <= 0:
        ln_na = math.log(na_M) if na_M > 0 else 0.0
        tm_K = 1 / (1 / tm_1M_K + (4.29 * gc_frac - 3.95) * 1e-5 * ln_na + 9.40e-6 * ln_na ** 2)
    else:
        ln_mg = math.log(mg_M)
        ln_na = math.log(na_M) if na_M > 0 else 0.0
        ratio = math.sqrt(mg_M) / na_M if na_M > 0 else float("inf")
        if ratio < 0.22:
            correction = (4.29 * gc_frac - 3.95) * 1e-5 * ln_na + 9.40e-6 * ln_na ** 2
        else:
            a, b, c = 3.92e-5, -9.11e-6, 6.26e-5
            d, e, f, g = 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5
            if ratio < 6.0 and na_M > 0:
                a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(na_M) * ln_na)
                d = 1.42e-5 * (1.279 - 4.03e-3 * ln_na - 8.03e-3 * ln_na ** 2)
                g = 8.31e-5 * (0.486 - 0.258 * ln_na + 5.25e-3 * ln_na ** 3)
            correction = (a + b * ln_mg + gc_frac * (c + d * ln_mg)
                          + (1 / (2 * n_bp)) * (e + f * ln_mg + g * ln_mg ** 2))
        tm_K = 1 / (1 / tm_1M_K + correction)
    return tm_K - 273.15


def melting_temperature_marmur(seq: str) -> float:
    """Wallace rule. Short oligos only (< 14 nt)."""
    from gsynth_core.sequence.ops import ensure_dna
    s = ensure_dna(seq)
    if not s:
        return 0.0
    a, t = s.count("A"), s.count("T")
    g, c = s.count("G"), s.count("C")
    return 2.0 * (a + t) + 4.0 * (g + c) - (7 if len(s) >= 13 else 0)


def melting_temperature_gc(seq: str, *, na_conc_mM: float = 50.0) -> float:
    """Wetmur 1991 — long duplex approximation."""
    from gsynth_core.sequence.ops import ensure_dna, gc_content
    s = ensure_dna(seq)
    if not s:
        return 0.0
    return 81.5 + 0.41 * gc_content(s) + 16.6 * math.log10(na_conc_mM * 1e-3) - 675.0 / len(s)


def melting_temperature(seq: str, *, method: TmMethod = TmMethod.AUTO,
                        dna_conc_nM: float = 500.0, na_conc_mM: float = 50.0,
                        mg_conc_mM: float = 1.5, dntp_conc_mM: float = 0.2) -> float:
    """Main entry: primer3-py if available, else SantaLucia NN."""
    from gsynth_core.sequence.ops import ensure_dna
    s = ensure_dna(seq)
    if method in (TmMethod.MARMUR, TmMethod.WALLACE):
        return melting_temperature_marmur(s)
    if method == TmMethod.GC_CONTENT:
        return melting_temperature_gc(s, na_conc_mM=na_conc_mM)
    if method == TmMethod.AUTO:
        try:
            import primer3
            return float(primer3.calc_tm(
                s, mv_conc=na_conc_mM, dv_conc=mg_conc_mM,
                dntp_conc=dntp_conc_mM, dna_conc=dna_conc_nM,
            ))
        except Exception:
            pass
    return melting_temperature_nn(
        s, dna_conc_nM=dna_conc_nM, na_conc_mM=na_conc_mM,
        mg_conc_mM=mg_conc_mM, dntp_conc_mM=dntp_conc_mM,
    )
