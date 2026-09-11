from __future__ import annotations

from typing import Final

from gsynth_engine.enzyme_table import ENZYME_TABLE

HIS_TAG: Final[str] = "CACCACCACCACCACCAC"
LEFT_LINKER: Final[str] = "GGTTCTTCT"
RIGHT_LINKER: Final[str] = "TCTTCTGGT"

CLEAVAGE_SITES: Final[dict[str, str]] = {
    "Thrombin":     "CTGGTGCCGCGTGGTTCT",
    "TEV":          "GAAAACCTGTATTTTCAGGGC",


    "Factor Xa":    "ATCGAAGGTCGT",
    "PreScission":  "CTGGAAGTGCTGTTCCAGGGCCCA",
    "Enterokinase": "GATGACGATGACAAG",
    "SUMO":         "CTGCAGGACTCAGAGG",
    "HRV 3C":       "CTGGAAGTTCTGTTCCAGGGGCCC",
}


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


ALL_ENZYMES: Final[dict[str, dict[str, object]]] = {
    **{name: dict(spec) for name, spec in ENZYME_TABLE.items()},
    **{
        name: {**ENZYME_TABLE.get(name, {}), **spec}
        for name, spec in RESTRICTION_ENZYMES.items()
    },
}


def _rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def left_remainders(enzyme: str) -> tuple[str, str]:

    info = ALL_ENZYMES[enzyme]
    site: str = info["recognition"]           # type: ignore[assignment]
    top: int = info["cut_top"]                # type: ignore[assignment]
    bottom: int = info["cut_bottom"]          # type: ignore[assignment]
    return site[top:], _rc(site[bottom:])


def supplies_start_codon(enzyme: str) -> bool:

    forward, _reverse = left_remainders(enzyme)
    return forward.endswith("ATG")


def right_remainders(enzyme: str) -> tuple[str, str]:

    info = ALL_ENZYMES[enzyme]
    site: str = info["recognition"]           # type: ignore[assignment]
    top: int = info["cut_top"]                # type: ignore[assignment]
    bottom: int = info["cut_bottom"]          # type: ignore[assignment]
    return site[:top], _rc(site[:bottom])


def overhang(enzyme: str) -> tuple[str, str]:

    info = ALL_ENZYMES[enzyme]
    site: str = info["recognition"]           # type: ignore[assignment]
    top: int = info["cut_top"]                # type: ignore[assignment]
    bottom: int = info["cut_bottom"]          # type: ignore[assignment]
    if top == bottom:
        return "", "blunt"
    lo, hi = min(top, bottom), max(top, bottom)
    return site[lo:hi], ("5'" if top < bottom else "3'")


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
