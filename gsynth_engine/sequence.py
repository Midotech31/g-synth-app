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


def melting_temperature(sequence: str, *, na_mM: float = 50.0) -> float:
    """Tm in °C.

    Wallace rule below 14 nt, where nearest-neighbour models are unreliable;
    the salt-adjusted GC formula above it. Good enough to compare oligos in
    the same design — order-critical work should use a dedicated tool.
    """
    seq = clean_dna(sequence)
    if not seq:
        return 0.0
    a, t = seq.count("A"), seq.count("T")
    g, c = seq.count("G"), seq.count("C")
    if len(seq) < 14:
        return 2.0 * (a + t) + 4.0 * (g + c)
    return (
        81.5
        + 16.6 * math.log10(max(na_mM, 1e-6) / 1000.0)
        + 0.41 * gc_content(seq)
        - 675.0 / len(seq)
    )
