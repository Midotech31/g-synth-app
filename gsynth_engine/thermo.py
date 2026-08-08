"""Melting temperatures by the nearest-neighbour model.

SnapGene, Geneious and Primer3 all compute Tm from stacking thermodynamics
rather than from base composition, because a sequence's Tm depends on the
order of its bases, not just how many G and C it contains. `AAGGCC` and
`AGAGCC` have identical GC content and different melting behaviour; a
composition formula cannot tell them apart.

Implemented here:

* **SantaLucia 1998** unified nearest-neighbour parameters (ΔH kcal/mol,
  ΔS cal/mol/K) with initiation terms for each end.
* **Owczarzy 2004/2008** salt correction, including the divalent-cation
  path used when Mg²⁺ is present — the case that matters for annealing
  buffers, where the monovalent formula overestimates Tm.
* The symmetry penalty for self-complementary duplexes.

References
    SantaLucia J. (1998) PNAS 95:1460–1465.
    Owczarzy R. et al. (2004) Biochemistry 43:3537–3554.
    Owczarzy R. et al. (2008) Biochemistry 47:5336–5353.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

from gsynth_engine.sequence import (
    clean_dna,
    gc_content,
    is_palindrome,
    reverse_complement,
)


@dataclass(frozen=True)
class BufferConditions:
    """The reaction a melting temperature refers to.

    A Tm quoted without its conditions is not actionable: the same oligo
    melts several degrees apart in a primer dilution and in the 25 µM
    annealing reaction this protocol prescribes. Every Tm the engine
    reports carries the conditions it was computed under.
    """

    name: str
    oligo_nM: float   #: total strand concentration Ct, not per-strand
    na_mM: float
    mg_mM: float = 0.0
    dntp_mM: float = 0.0

    @property
    def summary(self) -> str:
        """One line, for an order sheet footnote or a UI tooltip."""
        parts = [f"{self.oligo_nM / 1000:g} µM total strand", f"{self.na_mM:g} mM Na⁺"]
        if self.mg_mM:
            parts.append(f"{self.mg_mM:g} mM Mg²⁺")
        if self.dntp_mM:
            parts.append(f"{self.dntp_mM:g} mM dNTP")
        return ", ".join(parts)


#: The annealing reaction in `protocol.py`: 5 µL + 5 µL of 100 µM oligos in
#: 20 µL, so 25 µM of each strand and Ct = 50 µM, in 1× annealing buffer
#: (50 mM NaCl, no Mg²⁺). This is what the oligos actually experience, so it
#: is what the design pages report.
ANNEALING: Final[BufferConditions] = BufferConditions(
    name="annealing reaction", oligo_nM=50_000.0, na_mM=50.0
)

#: A conventional primer Tm, for comparison against a supplier's datasheet.
PRIMER: Final[BufferConditions] = BufferConditions(
    name="primer", oligo_nM=500.0, na_mM=50.0
)

#: SantaLucia 1998 unified parameters, ΔH (kcal/mol) and ΔS (cal/mol/K).
#: Keyed "top/bottom", both written 5'→3'. The ten entries below cover all
#: sixteen dinucleotides once dyad symmetry is applied in `_stack`.
_NN: Final[dict[str, tuple[float, float]]] = {
    "AA/TT": (-7.9, -22.2),
    "AT/AT": (-7.2, -20.4),
    "TA/TA": (-7.2, -21.3),
    "CA/TG": (-8.5, -22.7),
    "GT/AC": (-8.4, -22.4),
    "CT/AG": (-7.8, -21.0),
    "GA/TC": (-8.2, -22.2),
    "CG/CG": (-10.6, -27.2),
    "GC/GC": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}

#: Helix initiation, applied once per end according to that end's base.
_INIT_GC: Final[tuple[float, float]] = (0.1, -2.8)
_INIT_AT: Final[tuple[float, float]] = (2.3, 4.1)

#: Entropy penalty for a self-complementary duplex.
_SYMMETRY_DS: Final[float] = -1.4

_R: Final[float] = 1.9872  # gas constant, cal/(mol·K)


def _stack(dimer: str) -> tuple[float, float]:
    """Nearest-neighbour parameters for one dinucleotide step."""
    top = dimer.upper()
    bottom = reverse_complement(top)
    for key in (f"{top}/{bottom}", f"{bottom}/{top}"):
        if key in _NN:
            return _NN[key]
    raise KeyError(f"No nearest-neighbour parameters for {dimer!r}")


def duplex_thermodynamics(sequence: str) -> tuple[float, float]:
    """Total ΔH (kcal/mol) and ΔS (cal/mol/K) for a perfectly matched duplex."""
    seq = clean_dna(sequence)
    if len(seq) < 2:
        return 0.0, 0.0

    enthalpy = entropy = 0.0
    for i in range(len(seq) - 1):
        step_h, step_s = _stack(seq[i : i + 2])
        enthalpy += step_h
        entropy += step_s

    for end in (seq[0], seq[-1]):
        init_h, init_s = _INIT_GC if end in "GC" else _INIT_AT
        enthalpy += init_h
        entropy += init_s

    if is_palindrome(seq):
        entropy += _SYMMETRY_DS

    return enthalpy, entropy


def _salt_corrected(
    tm_kelvin: float, *, gc_fraction: float, n_phosphates: int,
    na_mM: float, mg_mM: float, dntp_mM: float,
) -> float:
    """Owczarzy salt correction, monovalent or divalent as appropriate."""
    # dNTPs chelate Mg²⁺; only the free fraction affects stability.
    free_mg_M = max(0.0, (mg_mM - dntp_mM)) * 1e-3
    na_M = max(na_mM, 1e-6) * 1e-3
    ln_na = math.log(na_M)

    monovalent = (4.29 * gc_fraction - 3.95) * 1e-5 * ln_na + 9.40e-6 * ln_na**2

    if free_mg_M <= 0:
        return 1.0 / (1.0 / tm_kelvin + monovalent)

    ratio = math.sqrt(free_mg_M) / na_M if na_M > 0 else float("inf")
    if ratio < 0.22:
        # Monovalent ions dominate.
        return 1.0 / (1.0 / tm_kelvin + monovalent)

    ln_mg = math.log(free_mg_M)
    a, b, c = 3.92e-5, -9.11e-6, 6.26e-5
    d, e, f, g = 1.42e-5, -4.82e-4, 5.25e-4, 8.31e-5
    if ratio < 6.0 and na_M > 0:
        # Both ion types contribute; Owczarzy 2008 eq. 19 adjusts a, d and g.
        a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(na_M) * ln_na)
        d = 1.42e-5 * (1.279 - 4.03e-3 * ln_na - 8.03e-3 * ln_na**2)
        g = 8.31e-5 * (0.486 - 0.258 * ln_na + 5.25e-3 * ln_na**3)

    divalent = (
        a
        + b * ln_mg
        + gc_fraction * (c + d * ln_mg)
        + (1.0 / (2.0 * max(n_phosphates, 1))) * (e + f * ln_mg + g * ln_mg**2)
    )
    return 1.0 / (1.0 / tm_kelvin + divalent)


def melting_temperature(
    sequence: str,
    *,
    conditions: BufferConditions = PRIMER,
    oligo_nM: float | None = None,
    na_mM: float | None = None,
    mg_mM: float | None = None,
    dntp_mM: float | None = None,
) -> float:
    """Tm in °C for a perfectly matched duplex.

    `conditions` sets the reaction; the individual keywords override single
    values within it. The default is a conventional primer Tm — pass
    `conditions=ANNEALING` for the reaction the G-Synth protocol prescribes.

    Sequences shorter than 8 nt fall back to the Wallace rule, where the
    nearest-neighbour model is not parameterised — junction overhangs of
    4–8 nt land in that range deliberately.
    """
    oligo_nM = conditions.oligo_nM if oligo_nM is None else oligo_nM
    na_mM = conditions.na_mM if na_mM is None else na_mM
    mg_mM = conditions.mg_mM if mg_mM is None else mg_mM
    dntp_mM = conditions.dntp_mM if dntp_mM is None else dntp_mM

    seq = clean_dna(sequence)
    if not seq:
        return 0.0
    if len(seq) < 8:
        a_t = seq.count("A") + seq.count("T")
        g_c = seq.count("G") + seq.count("C")
        return 2.0 * a_t + 4.0 * g_c

    enthalpy, entropy = duplex_thermodynamics(seq)

    # Self-complementary strands are their own partner, so the concentration
    # term uses Ct rather than Ct/4.
    total_strand_M = oligo_nM * 1e-9
    divisor = 1.0 if is_palindrome(seq) else 4.0

    tm_kelvin = (enthalpy * 1000.0) / (
        entropy + _R * math.log(total_strand_M / divisor)
    )
    tm_kelvin = _salt_corrected(
        tm_kelvin,
        gc_fraction=gc_content(seq) / 100.0,
        n_phosphates=len(seq) - 1,
        na_mM=na_mM,
        mg_mM=mg_mM,
        dntp_mM=dntp_mM,
    )
    return tm_kelvin - 273.15
