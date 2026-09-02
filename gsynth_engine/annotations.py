"""Conservative recognition of common sequence features.

This module deliberately reports *exact motif matches*, not inferred biological
function.  A matching T7 promoter motif is useful evidence for annotation, but
it does not prove that the promoter is active in a particular construct or
host.  Callers therefore present these results for review before saving them.

Coordinates follow the engine convention: zero-based, half-open.  For a
circular record, ``end`` may be greater than the sequence length when a motif
crosses the origin.
"""
from __future__ import annotations

from dataclasses import dataclass

from gsynth_engine.constants import (
    CLEAVAGE_SITES,
    HIS_TAG,
    LEFT_LINKER,
    RIGHT_LINKER,
)
from gsynth_engine.sequence import reverse_complement, validate_dna


@dataclass(frozen=True)
class CommonMotif:
    """One curated exact DNA spelling and its proposed map annotation."""

    name: str
    sequence: str
    feature_type: str
    color: str
    basis: str


COMMON_MOTIFS: tuple[CommonMotif, ...] = (
    CommonMotif(
        "T7 promoter",
        "TAATACGACTCACTATAGGG",
        "promoter",
        "#C97634",
        "Exact match to the canonical T7 promoter motif.",
    ),
    CommonMotif(
        "lac operator",
        "AATTGTGAGCGGATAACAATT",
        "protein_bind",
        "#4A7C59",
        "Exact match to the canonical lac operator core.",
    ),
    CommonMotif(
        "T7 terminator",
        "CAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAG",
        "terminator",
        "#9E3D3D",
        "Exact match to the T7 terminator carried by the bundled pET record.",
    ),
    CommonMotif(
        "6×His tag",
        HIS_TAG,
        "misc_feature",
        "#0E6E77",
        "Exact match to the G-Synth six-histidine cassette DNA spelling.",
    ),
    CommonMotif(
        "Gly-Ser-Ser linker",
        LEFT_LINKER,
        "misc_feature",
        "#78889B",
        "Exact match to the G-Synth GSS linker DNA spelling.",
    ),
    CommonMotif(
        "Ser-Ser-Gly linker",
        RIGHT_LINKER,
        "misc_feature",
        "#78889B",
        "Exact match to the G-Synth SSG linker DNA spelling.",
    ),
    *tuple(
        CommonMotif(
            f"{name} site",
            sequence,
            "misc_feature",
            "#6A4C93",
            f"Exact match to the G-Synth {name} recognition-site DNA spelling.",
        )
        for name, sequence in CLEAVAGE_SITES.items()
    ),
)


def _occurrences(haystack: str, needle: str) -> list[int]:
    """Return overlapping occurrences, including tandem repeated motifs."""
    positions: list[int] = []
    start = 0
    while True:
        found = haystack.find(needle, start)
        if found < 0:
            return positions
        positions.append(found)
        start = found + 1


def _same_known_feature(existing: dict, motif: CommonMotif, start: int, end: int) -> bool:
    """Suppress a proposal already covered by the same named feature.

    Imported records often include two flanking bases around a canonical
    core.  Requiring identical coordinates would therefore propose a duplicate
    lac operator inside an annotation that already says ``lac operator``.
    """
    existing_name = str(existing.get("name", "")).casefold().replace("×", "x")
    motif_name = motif.name.casefold().replace("×", "x")
    same_name = motif_name in existing_name or existing_name in motif_name
    try:
        existing_start = int(existing.get("start", -1))
        existing_end = int(existing.get("end", -1))
    except (TypeError, ValueError):
        return False
    overlap = max(0, min(existing_end, end) - max(existing_start, start))
    # Imported feature boundaries commonly include or omit one flanking base.
    # A same-named feature covering at least 80% of the exact motif is already
    # the user's annotation; proposing a second bar would add noise, not truth.
    substantially_same_span = overlap / max(1, end - start) >= 0.8
    return same_name and substantially_same_span


def detect_common_features(
    sequence: str,
    *,
    circular: bool = False,
    existing: list[dict] | None = None,
) -> list[dict]:
    """Return reviewable exact matches to the curated common-motif library.

    Both strands are searched.  Palindromic motifs are emitted once, and a
    circular search reports a wrapped coordinate only once at its true start.
    Existing features with the same name that cover the match are omitted.
    """
    dna = validate_dna(sequence)
    length = len(dna)
    existing = existing or []
    matches: list[dict] = []
    seen: set[tuple[str, int, int]] = set()

    for motif in COMMON_MOTIFS:
        motif_dna = validate_dna(motif.sequence, field=f"{motif.name} motif")
        extension = dna[: len(motif_dna) - 1] if circular else ""
        searchable = dna + extension
        strands = ((1, motif_dna), (-1, reverse_complement(motif_dna)))
        for direction, query in strands:
            for start in _occurrences(searchable, query):
                if start >= length:
                    continue
                end = start + len(query)
                if not circular and end > length:
                    continue
                key = (motif.name, start, end)
                if key in seen:
                    continue
                seen.add(key)
                if any(_same_known_feature(item, motif, start, end) for item in existing):
                    continue
                matches.append({
                    "annotation": {
                        "name": motif.name,
                        "type": motif.feature_type,
                        "start": start,
                        "end": end,
                        "direction": direction,
                        "color": motif.color,
                    },
                    "matched_sequence": query,
                    "basis": motif.basis,
                })

    return sorted(
        matches,
        key=lambda match: (
            match["annotation"]["start"],
            match["annotation"]["end"],
            match["annotation"]["name"],
        ),
    )
