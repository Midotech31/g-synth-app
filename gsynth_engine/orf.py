"""Translation and six-frame ORF analysis for linear DNA.

Coordinates are always reported on the sequence the user supplied: 0-based,
half-open, and left-to-right even when the ORF is on the reverse strand.  A
negative frame states the strand and the phase together (``-1`` through
``-3``), matching the conventional six-frame view used at the bench.
"""
from __future__ import annotations

from dataclasses import dataclass

from gsynth_engine.cloning import translate
from gsynth_engine.constants import STOP_CODONS
from gsynth_engine.sequence import gc_content, reverse_complement, validate_dna


@dataclass(frozen=True)
class TranslationFrame:
    frame: int
    strand: str
    offset: int
    dna: str
    protein: str
    first_atg: int | None
    protein_from_first_atg: str


@dataclass(frozen=True)
class ORF:
    index: int
    frame: int
    strand: str
    start: int
    end: int
    length_nt: int
    amino_acids: int
    dna: str
    protein: str
    start_codon: str
    stop_codon: str


@dataclass(frozen=True)
class SequenceAnalysis:
    sequence: str
    reverse_complement: str
    length: int
    gc: float
    frames: tuple[TranslationFrame, ...]
    orfs: tuple[ORF, ...]
    minimum_codons: int


def _frame_payload(sequence: str, *, reverse: bool, offset: int) -> TranslationFrame:
    strand_sequence = reverse_complement(sequence) if reverse else sequence
    framed = strand_sequence[offset:]
    first = next(
        (position for position in range(0, len(framed) - 2, 3)
         if framed[position:position + 3] == "ATG"),
        None,
    )
    if first is None:
        first_atg = None
        from_first = ""
    else:
        strand_at = offset + first
        first_atg = len(sequence) - strand_at - 3 if reverse else strand_at
        from_first = translate(framed[first:])

    frame = -(offset + 1) if reverse else offset + 1
    return TranslationFrame(
        frame=frame,
        strand="reverse" if reverse else "forward",
        offset=offset,
        dna=framed,
        protein=translate(framed),
        first_atg=first_atg,
        protein_from_first_atg=from_first,
    )


def _find_orfs_on_strand(
    sequence: str, *, reverse: bool, minimum_codons: int,
) -> list[ORF]:
    strand_sequence = reverse_complement(sequence) if reverse else sequence
    length = len(sequence)
    found: list[ORF] = []

    for offset in range(3):
        starts: list[int] = []
        for position in range(offset, len(strand_sequence) - 2, 3):
            codon = strand_sequence[position:position + 3]
            if codon == "ATG":
                starts.append(position)
                continue
            if codon not in STOP_CODONS or not starts:
                continue

            stop = position + 3
            for start in starts:
                coding_codons = (stop - start) // 3 - 1
                if coding_codons < minimum_codons:
                    continue
                dna = strand_sequence[start:stop]
                if reverse:
                    top_start, top_end = length - stop, length - start
                else:
                    top_start, top_end = start, stop
                found.append(ORF(
                    index=0,
                    frame=-(offset + 1) if reverse else offset + 1,
                    strand="reverse" if reverse else "forward",
                    start=top_start,
                    end=top_end,
                    length_nt=len(dna),
                    amino_acids=coding_codons,
                    dna=dna,
                    protein=translate(dna)[:-1],
                    start_codon="ATG",
                    stop_codon=codon,
                ))
            starts = []

    found.sort(key=lambda item: (-item.amino_acids, item.start, item.frame))
    return found


def analyse_sequence(sequence: str, *, minimum_codons: int = 10) -> SequenceAnalysis:
    """Translate all six frames and return every complete ATG-to-stop ORF.

    ``minimum_codons`` counts translated amino acids and excludes the terminal
    stop. Nested start codons are retained as distinct candidate ORFs because
    an internal methionine can be the biologically used start.
    """
    if not 1 <= minimum_codons <= 100_000:
        raise ValueError("minimum_codons must be between 1 and 100000")
    seq = validate_dna(sequence, field="DNA sequence")
    reverse = reverse_complement(seq)
    frames = tuple(
        _frame_payload(seq, reverse=is_reverse, offset=offset)
        for is_reverse in (False, True)
        for offset in range(3)
    )
    records = (
        _find_orfs_on_strand(seq, reverse=False, minimum_codons=minimum_codons)
        + _find_orfs_on_strand(seq, reverse=True, minimum_codons=minimum_codons)
    )
    records.sort(key=lambda item: (-item.amino_acids, item.start, item.frame))
    indexed = tuple(
        ORF(**{**record.__dict__, "index": index})
        for index, record in enumerate(records, start=1)
    )
    return SequenceAnalysis(
        sequence=seq,
        reverse_complement=reverse,
        length=len(seq),
        gc=round(gc_content(seq), 1),
        frames=frames,
        orfs=indexed,
        minimum_codons=minimum_codons,
    )
