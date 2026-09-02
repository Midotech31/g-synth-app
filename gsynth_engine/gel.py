"""In-silico agarose-gel inputs derived from sequence and restriction cuts.

The engine reports fragment sizes.  It does not fabricate experimental lane
brightness, background, smearing or topology-dependent migration.  The UI can
plot these sizes as a clearly labelled prediction and choose a reference
ladder that brackets the expected bands.
"""
from __future__ import annotations

from typing import Final

from gsynth_engine.cloning import _cut_positions, find_sites
from gsynth_engine.sequence import SequenceError, validate_dna

GEL_LADDERS: Final[dict[str, dict[str, object]]] = {
    "100-bp": {
        "name": "100 bp reference ladder",
        "bands": (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500),
        "range": "100–1,500 bp",
    },
    "1-kb": {
        "name": "1 kb reference ladder",
        "bands": (500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000),
        "range": "500–10,000 bp",
    },
    "broad-range": {
        "name": "Broad-range DNA ladder",
        "bands": (100, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000, 4000, 5000, 7000, 10000),
        "range": "100–10,000 bp",
    },
}


def recommended_ladder(fragment_sizes: list[int]) -> str:
    """Choose the narrowest generic ladder that brackets all useful bands."""
    if not fragment_sizes:
        return "broad-range"
    maximum = max(fragment_sizes)
    minimum = min(fragment_sizes)
    if minimum >= 100 and maximum <= 1500:
        return "100-bp"
    if minimum >= 500 and maximum <= 10000:
        return "1-kb"
    return "broad-range"


def restriction_digest_sizes(
    sequence: str,
    enzymes: list[str] | tuple[str, ...],
    *,
    circular: bool = True,
) -> list[int]:
    """Return all fragment lengths from a complete in-silico digest.

    Cut coordinates use each enzyme's top-strand cut offset.  Duplicate cut
    coordinates are collapsed, as two enzymes cutting the same phosphodiester
    bond do not create a zero-length gel fragment.
    """
    dna = validate_dna(sequence)
    if not enzymes:
        raise SequenceError("Choose at least one restriction enzyme for a digest simulation.")

    cuts: set[int] = set()
    for enzyme in enzymes:
        positions = find_sites(dna, enzyme, circular=circular)
        if not positions:
            raise SequenceError(f"{enzyme} does not cut this sequence.")
        cuts.update(_cut_positions(enzyme, position)[0] % len(dna) for position in positions)

    ordered = sorted(cuts)
    if circular:
        if len(ordered) == 1:
            return [len(dna)]
        sizes = [right - left for left, right in zip(ordered, ordered[1:], strict=False)]
        sizes.append(len(dna) - ordered[-1] + ordered[0])
    else:
        boundaries = [0, *ordered, len(dna)]
        sizes = [right - left for left, right in zip(boundaries, boundaries[1:], strict=False)]

    return sorted((size for size in sizes if size > 0), reverse=True)


def ladder_payload() -> list[dict]:
    """JSON-ready ladder catalogue with every marker size visible."""
    return [
        {"key": key, "name": entry["name"], "bands": list(entry["bands"]), "range": entry["range"]}
        for key, entry in GEL_LADDERS.items()
    ]
