"""Biological constants used by the G-Synth design engine.

These values are the G-Synth standard cassette. They are reproduced exactly
as specified — the golden tests depend on them, and changing one changes the
oligos the lab orders.
"""
from __future__ import annotations

from typing import Final

from gsynth_engine.enzyme_table import ENZYME_TABLE

# ── Standard purification / cleavage cassette ────────────────────────────────
HIS_TAG: Final[str] = "CACCACCACCACCACCAC"          # 6× His
LEFT_LINKER: Final[str] = "GGTTCTTCT"               # flexible Gly-Ser-Ser
RIGHT_LINKER: Final[str] = "TCTTCTGGT"              # flexible Ser-Ser-Gly

CLEAVAGE_SITES: Final[dict[str, str]] = {
    "Thrombin":     "CTGGTGCCGCGTGGTTCT",           # LVPR/GS
    "TEV":          "GAAAACCTGTATTTTCAGGGC",        # ENLYFQ/G
    "Factor Xa":    "ATCGAAGGTCGT",                 # IEGR/
    "PreScission":  "CTGGAAGTGCTGTTCCAGGGCCCA",     # LEVLFQ/GP
    "Enterokinase": "GATGACGATGACAAG",              # DDDDK/
    "SUMO":         "CTGCAGGACTCAGAGG",
    "HRV 3C":       "CTGGAAGTTCTGTTCCAGGGGCCC",     # LEVLFQ/GP
}

# ── Restriction enzymes ──────────────────────────────────────────────────────
# One source of truth per enzyme: the recognition site and where the two
# strands are cut inside it. `cut_top` / `cut_bottom` are offsets from the
# start of the site, in top-strand coordinates.
#
#     NdeI   CA^TATG   cut_top = 2, cut_bottom = 4  → 5'-TA overhang
#     XhoI   C^TCGAG   cut_top = 1, cut_bottom = 5  → 5'-TCGA overhang
#     KpnI   GGTAC^C   cut_top = 5, cut_bottom = 1  → 3'-GTAC overhang
#     SmaI   CCC^GGG   cut_top = cut_bottom = 3     → blunt
#
# What each oligo must carry is *derived* from these (see `left_remainders`
# and `right_remainders`) rather than stored, because it differs depending on
# whether the enzyme sits at the left or the right end of the insert. Storing
# a single pair per enzyme — as G-Synth 2.x did — silently produced
# mismatched duplexes for every enzyme except the validated NdeI/XhoI pair.
#
# NdeI is the important special case: its site CATATG *contains* the ATG
# start codon, so the forward oligo carries "TATG" — the T completes CATATG
# once ligated, and the A of the 5'-TA overhang is also the A of the ATG.
# A separate ATG must therefore never be added for NdeI.
RESTRICTION_ENZYMES: Final[dict[str, dict[str, object]]] = {
    "NdeI":    {"recognition": "CATATG",   "cut_top": 2, "cut_bottom": 4},
    "XhoI":    {"recognition": "CTCGAG",   "cut_top": 1, "cut_bottom": 5},
    "EcoRI":   {"recognition": "GAATTC",   "cut_top": 1, "cut_bottom": 5},
    "BamHI":   {"recognition": "GGATCC",   "cut_top": 1, "cut_bottom": 5},
    "HindIII": {"recognition": "AAGCTT",   "cut_top": 1, "cut_bottom": 5},
    "SalI":    {"recognition": "GTCGAC",   "cut_top": 1, "cut_bottom": 5},
    "XbaI":    {"recognition": "TCTAGA",   "cut_top": 1, "cut_bottom": 5},
    "NcoI":    {"recognition": "CCATGG",   "cut_top": 1, "cut_bottom": 5},
    "BglII":   {"recognition": "AGATCT",   "cut_top": 1, "cut_bottom": 5},
    "SpeI":    {"recognition": "ACTAGT",   "cut_top": 1, "cut_bottom": 5},
    "MluI":    {"recognition": "ACGCGT",   "cut_top": 1, "cut_bottom": 5},
    "NotI":    {"recognition": "GCGGCCGC", "cut_top": 2, "cut_bottom": 6},
    "KpnI":    {"recognition": "GGTACC",   "cut_top": 5, "cut_bottom": 1},
    "SacI":    {"recognition": "GAGCTC",   "cut_top": 5, "cut_bottom": 1},
    "PstI":    {"recognition": "CTGCAG",   "cut_top": 5, "cut_bottom": 1},
    "ApaI":    {"recognition": "GGGCCC",   "cut_top": 5, "cut_bottom": 1},
    "SmaI":    {"recognition": "CCCGGG",   "cut_top": 3, "cut_bottom": 3},
    "EcoRV":   {"recognition": "GATATC",   "cut_top": 3, "cut_bottom": 3},
    "SspI":    {"recognition": "AATATT",   "cut_top": 3, "cut_bottom": 3},
}


#: Every enzyme G-Synth can reason about: the curated set above, plus the
#: wide REBASE table. Use this to ask "what cuts here". `RESTRICTION_ENZYMES`
#: remains the list a user *picks a cloning pair from* — a dropdown of a
#: hundred names is worse than nineteen, and this lab uses nineteen.
ALL_ENZYMES: Final[dict[str, dict[str, object]]] = {
    **{name: dict(spec) for name, spec in ENZYME_TABLE.items()},
    **RESTRICTION_ENZYMES,          # the curated entries win on name collision
}


def _rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def left_remainders(enzyme: str) -> tuple[str, str]:
    """What each oligo carries when this enzyme is at the LEFT end.

    Returns (forward_prefix, reverse_suffix): the forward oligo starts with
    the first, the reverse oligo ends with the second. The insert is the
    right-hand piece of the cut, so it keeps everything after the cut.
    """
    info = ALL_ENZYMES[enzyme]
    site: str = info["recognition"]           # type: ignore[assignment]
    top: int = info["cut_top"]                # type: ignore[assignment]
    bottom: int = info["cut_bottom"]          # type: ignore[assignment]
    return site[top:], _rc(site[bottom:])


def right_remainders(enzyme: str) -> tuple[str, str]:
    """What each oligo carries when this enzyme is at the RIGHT end.

    Returns (forward_suffix, reverse_prefix). The insert is the left-hand
    piece of the cut, so it keeps everything before the cut.
    """
    info = ALL_ENZYMES[enzyme]
    site: str = info["recognition"]           # type: ignore[assignment]
    top: int = info["cut_top"]                # type: ignore[assignment]
    bottom: int = info["cut_bottom"]          # type: ignore[assignment]
    return site[:top], _rc(site[:bottom])


def overhang(enzyme: str) -> tuple[str, str]:
    """The sticky end this enzyme leaves: (sequence, "5'" | "3'" | "blunt").

    Sequence is given in top-strand sense.
    """
    info = ALL_ENZYMES[enzyme]
    site: str = info["recognition"]           # type: ignore[assignment]
    top: int = info["cut_top"]                # type: ignore[assignment]
    bottom: int = info["cut_bottom"]          # type: ignore[assignment]
    if top == bottom:
        return "", "blunt"
    lo, hi = min(top, bottom), max(top, bottom)
    return site[lo:hi], ("5'" if top < bottom else "3'")

#: Pairs routinely used for pET cloning, offered first in the UI.
COMMON_ENZYME_PAIRS: Final[tuple[str, ...]] = (
    "NdeI / XhoI",
    "NdeI / EcoRI",
    "NcoI / XhoI",
    "BamHI / EcoRI",
    "BamHI / XhoI",
    "EcoRI / SalI",
    "NdeI / HindIII",
)

STOP_CODONS: Final[frozenset[str]] = frozenset({"TAA", "TAG", "TGA"})
