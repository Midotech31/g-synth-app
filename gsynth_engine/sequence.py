from __future__ import annotations

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


class SequenceError(ValueError):
    pass


def clean_dna(sequence: str) -> str:

    lines = [ln for ln in sequence.splitlines() if not ln.lstrip().startswith(">")]
    return "".join(ch for ch in "".join(lines).upper() if ch.isalpha())


def validate_dna(sequence: str, *, field: str = "sequence") -> str:

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

    if not sequence:
        return 0.0
    gc = sum(1 for base in sequence.upper() if base in "GC")
    return 100.0 * gc / len(sequence)


def is_palindrome(sequence: str) -> bool:

    seq = sequence.upper()
    return bool(seq) and seq == reverse_complement(seq)


def longest_homopolymer(sequence: str) -> int:

    best = run = 0
    previous = ""
    for base in sequence.upper():
        run = run + 1 if base == previous else 1
        previous = base
        best = max(best, run)
    return best


def melting_temperature(sequence: str, **kwargs) -> float:

    from gsynth_engine.thermo import melting_temperature as _nn

    return _nn(sequence, **kwargs)
