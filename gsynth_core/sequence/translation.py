"""Translation and 6-frame ORF finder.

Fixes G-Synth 2.x bugs:
- 6-frame search (not just 3 forward frames)
- Reverse-strand coordinate conversion
- Optional alternative start codons GTG/TTG
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Literal

from gsynth_core.constants import (
    GENETIC_CODE, START_CODONS_PROKARYOTIC, START_CODONS_STANDARD, STOP_CODONS,
)
from gsynth_core.sequence.ops import clean_dna, ensure_dna, reverse_complement

Strand = Literal["+", "-"]


@dataclass(frozen=True, slots=True)
class ORF:
    start: int
    end: int
    strand: Strand
    frame: int
    protein: str
    start_codon: str
    stop_codon: str | None

    @property
    def length_nt(self) -> int:
        return self.end - self.start

    @property
    def length_aa(self) -> int:
        return len(self.protein)

    def __repr__(self) -> str:
        return (f"ORF({self.start}-{self.end} {self.strand}{self.frame}, "
                f"{self.length_aa} aa, start={self.start_codon}, stop={self.stop_codon})")


def translate(seq: str, *, frame: int = 0, to_stop: bool = False, unknown_aa: str = "X") -> str:
    """Translate DNA → protein."""
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0, 1 or 2, got {frame}")
    cleaned = clean_dna(seq, allow_ambiguous=True)
    out = []
    for i in range(frame, len(cleaned) - 2, 3):
        aa = GENETIC_CODE.get(cleaned[i:i+3], unknown_aa)
        if to_stop and aa == "*":
            break
        out.append(aa)
    return "".join(out)


def six_frame_translate(seq: str, *, to_stop: bool = False) -> dict[str, str]:
    """All 6 reading frames."""
    cleaned = ensure_dna(seq, allow_ambiguous=True)
    rc = reverse_complement(cleaned)
    return {f"+{i+1}": translate(cleaned, frame=i, to_stop=to_stop) for i in range(3)} | {
        f"-{i+1}": translate(rc, frame=i, to_stop=to_stop) for i in range(3)}


def find_orfs(seq: str, *, min_length_aa: int = 30, both_strands: bool = True,
              allow_alt_starts: bool = False, require_stop: bool = True) -> list[ORF]:
    """Find all ORFs across 6 frames with correct coordinate conversion."""
    cleaned = ensure_dna(seq, allow_ambiguous=True)
    starts = START_CODONS_PROKARYOTIC if allow_alt_starts else START_CODONS_STANDARD
    results: list[ORF] = []
    strands: list[tuple[Strand, str]] = [("+", cleaned)]
    if both_strands:
        strands.append(("-", reverse_complement(cleaned)))
    L = len(cleaned)
    for strand, working in strands:
        for frame in range(3):
            i = frame
            n = len(working)
            while i <= n - 3:
                codon = working[i:i+3]
                if codon in starts:
                    protein_chars = [GENETIC_CODE.get(codon, "X")]
                    stop_codon: str | None = None
                    j = i + 3
                    while j <= n - 3:
                        nxt = working[j:j+3]
                        if nxt in STOP_CODONS:
                            stop_codon = nxt
                            break
                        protein_chars.append(GENETIC_CODE.get(nxt, "X"))
                        j += 3
                    orf_end_local = j + 3 if stop_codon is not None else j
                    if stop_codon is None and require_stop:
                        break
                    if len(protein_chars) >= min_length_aa:
                        if strand == "+":
                            s_fwd, e_fwd = i, orf_end_local
                        else:
                            s_fwd, e_fwd = L - orf_end_local, L - i
                        results.append(ORF(
                            start=s_fwd, end=e_fwd, strand=strand, frame=frame+1,
                            protein="".join(protein_chars),
                            start_codon=codon, stop_codon=stop_codon,
                        ))
                    i = orf_end_local
                else:
                    i += 3
    results.sort(key=lambda o: (o.start, o.strand))
    return results


def reverse_translate(protein: str, *, organism: str = "e_coli",
                      method: Literal["most_frequent", "random", "weighted"] = "most_frequent",
                      seed: int | None = None) -> str:
    """Protein → DNA. For real optimization use gsynth_core.codon.optimize."""
    import random as _r
    from gsynth_core.constants import CODON_TABLES
    if organism not in CODON_TABLES:
        raise ValueError(f"Unknown organism {organism!r}")
    table = CODON_TABLES[organism]
    rng = _r.Random(seed)
    out: list[str] = []
    for aa in protein.upper():
        if aa not in table:
            if aa == "*":
                t = table["*"]
                out.append(max(t, key=t.__getitem__))
            else:
                out.append("NNN")
            continue
        codons = table[aa]
        if method == "most_frequent":
            out.append(max(codons, key=codons.__getitem__))
        elif method == "random":
            out.append(rng.choice(list(codons)))
        elif method == "weighted":
            ks = list(codons); ws = list(codons.values())
            out.append(rng.choices(ks, weights=ws, k=1)[0])
        else:
            raise ValueError(f"unknown method {method!r}")
    return "".join(out)
