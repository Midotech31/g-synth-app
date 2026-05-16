"""DNA sequence operations — IUPAC-aware, case-preserving."""
from __future__ import annotations
import re
from dataclasses import dataclass

from gsynth_core.constants import COMPLEMENT_TABLE, IUPAC_DNA_ALPHABET


@dataclass(frozen=True, slots=True)
class ValidationResult:
    ok: bool
    cleaned: str
    message: str | None = None
    removed_chars: int = 0


_STRICT_DNA_RE = re.compile(r"[^ACGTacgt]")
_IUPAC_DNA_RE = re.compile(f"[^{IUPAC_DNA_ALPHABET}{IUPAC_DNA_ALPHABET.lower()}]")
_WS_RE = re.compile(r"\s+")


def clean_dna(seq: str, *, allow_ambiguous: bool = False, preserve_case: bool = False) -> str:
    """Strip whitespace and non-base characters."""
    if not seq:
        return ""
    s = _WS_RE.sub("", seq)
    if not preserve_case:
        s = s.upper()
    return (_IUPAC_DNA_RE if allow_ambiguous else _STRICT_DNA_RE).sub("", s)


def validate_dna(seq: str, *, allow_empty: bool = False, allow_ambiguous: bool = False,
                 min_length: int = 0, max_length: int | None = None) -> ValidationResult:
    """Validate without raising."""
    if not seq:
        return (ValidationResult(True, "") if allow_empty
                else ValidationResult(False, "", "Sequence is empty"))
    trimmed = _WS_RE.sub("", seq)
    cleaned = clean_dna(seq, allow_ambiguous=allow_ambiguous)
    removed = len(trimmed) - len(cleaned)
    if not cleaned:
        return ValidationResult(False, "", "Sequence contains no valid DNA characters", removed)
    if len(cleaned) < min_length:
        return ValidationResult(False, cleaned, f"Too short ({len(cleaned)} < {min_length})", removed)
    if max_length is not None and len(cleaned) > max_length:
        return ValidationResult(False, cleaned, f"Too long ({len(cleaned)} > {max_length})", removed)
    msg = f"Removed {removed} invalid character{'s' if removed != 1 else ''}" if removed else None
    return ValidationResult(True, cleaned, msg, removed)


def ensure_dna(seq: str, **kw) -> str:
    """Strict variant that raises on failure."""
    r = validate_dna(seq, **kw)
    if not r.ok:
        raise ValueError(r.message or "invalid DNA")
    return r.cleaned


def is_dna(seq: str, *, allow_ambiguous: bool = False) -> bool:
    if not seq:
        return False
    alpha = IUPAC_DNA_ALPHABET if allow_ambiguous else "ACGT"
    return all(c.upper() in alpha for c in seq)


def complement(seq: str) -> str:
    """Watson-Crick complement supporting full IUPAC; preserves case."""
    return seq.translate(COMPLEMENT_TABLE)


def reverse_complement(seq: str) -> str:
    return complement(seq)[::-1]


def gc_content(seq: str, *, as_fraction: bool = False) -> float:
    """GC content. Ambiguous bases ignored."""
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    if gc + at == 0:
        return 0.0
    frac = gc / (gc + at)
    return frac if as_fraction else frac * 100.0


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError(f"length mismatch: {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b) if x.upper() != y.upper())


def is_palindrome(seq: str) -> bool:
    s = seq.upper()
    return s == reverse_complement(s)


def _base_set(c: str) -> set[str]:
    from gsynth_core.constants import IUPAC_DNA
    return set(IUPAC_DNA.get(c.upper(), frozenset()))


def iupac_match(pattern: str, target: str) -> bool:
    """True if every position of target is compatible with pattern."""
    if len(pattern) != len(target):
        return False
    for p, t in zip(pattern.upper(), target.upper()):
        ps, ts = _base_set(p), _base_set(t)
        if not ps or not ts or not (ps & ts):
            return False
    return True


def find_motif(seq: str, motif: str, *, allow_iupac: bool = True) -> list[int]:
    """Return all 0-based start positions of motif on forward strand."""
    if not motif or not seq:
        return []
    if not allow_iupac:
        return [m.start() for m in re.finditer(f"(?={re.escape(motif)})", seq.upper())]
    n, m = len(seq), len(motif)
    return [i for i in range(n - m + 1) if iupac_match(motif, seq[i:i + m])]


def find_motif_both_strands(seq: str, motif: str, *, allow_iupac: bool = True) -> dict[str, list[int]]:
    """Forward and reverse strand motif search."""
    fwd = find_motif(seq, motif, allow_iupac=allow_iupac)
    rc_motif = reverse_complement(motif)
    rev = [] if rc_motif == motif else find_motif(seq, rc_motif, allow_iupac=allow_iupac)
    return {"+": fwd, "-": rev}
