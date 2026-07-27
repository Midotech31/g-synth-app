"""Small sequence utilities. Deliberately dependency-free."""
from __future__ import annotations

import math

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


class SequenceError(ValueError):
    """Raised for input a user can fix, with a message written for them."""


def clean_dna(sequence: str) -> str:
    """Strip whitespace, digits and FASTA headers; upper-case the rest."""
    lines = [ln for ln in sequence.splitlines() if not ln.lstrip().startswith(">")]
    return "".join(ch for ch in "".join(lines).upper() if ch.isalpha())


def validate_dna(sequence: str, *, field: str = "sequence") -> str:
    """Return a clean A/C/G/T string or raise with the offending characters."""
    cleaned = clean_dna(sequence)
    if not cleaned:
        raise SequenceError(f"The {field} is empty.")
    invalid = sorted(set(cleaned) - set("ACGT"))
    if invalid:
        raise SequenceError(
            f"The {field} contains characters that are not A, C, G or T: "
            + ", ".join(invalid)
        )
    return cleaned


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)


def gc_content(sequence: str) -> float:
    """GC percentage, 0.0 for an empty sequence."""
    if not sequence:
        return 0.0
    gc = sum(1 for base in sequence.upper() if base in "GC")
    return 100.0 * gc / len(sequence)


def is_palindrome(sequence: str) -> bool:
    """True when a sequence is its own reverse complement.

    A palindromic overhang anneals to itself, so two copies of the same
    fragment can ligate to each other — the classic silent cause of a
    scrambled assembly.
    """
    seq = sequence.upper()
    return bool(seq) and seq == reverse_complement(seq)


def longest_homopolymer(sequence: str) -> int:
    """Length of the longest single-base run."""
    best = run = 0
    previous = ""
    for base in sequence.upper():
        run = run + 1 if base == previous else 1
        previous = base
        best = max(best, run)
    return best


def melting_temperature(sequence: str, **kwargs) -> float:
    """Tm in °C — delegates to the nearest-neighbour model.

    Kept here so existing callers keep working, but the calculation lives in
    `gsynth_engine.thermo`: composition formulas cannot distinguish sequences
    with the same GC content and different stacking, which is exactly what
    oligo design needs to see. Imported lazily to avoid a circular import.
    """
    from gsynth_engine.thermo import melting_temperature as _nn

    return _nn(sequence, **kwargs)
