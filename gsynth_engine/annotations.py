from __future__ import annotations

from dataclasses import dataclass

from gsynth_engine.constants import (
    CLEAVAGE_SITES,
    HIS_TAG,
    LEFT_LINKER,
    RIGHT_LINKER,
)
from gsynth_engine.sequence import SequenceError, clean_dna, reverse_complement, validate_dna


@dataclass(frozen=True)
class CommonMotif:


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
        "Shine-Dalgarno candidate",
        "AAGGAG",
        "RBS",
        "#B8860B",
        "Exact match to a commonly used bacterial Shine-Dalgarno motif; positional review is required.",
    ),
    CommonMotif(
        "Shine-Dalgarno candidate",
        "AGGAGG",
        "RBS",
        "#B8860B",
        "Exact match to a commonly used bacterial Shine-Dalgarno motif; positional review is required.",
    ),
    CommonMotif(
        "T7 terminator",
        reverse_complement("CAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAG"),
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

    positions: list[int] = []
    start = 0
    while True:
        found = haystack.find(needle, start)
        if found < 0:
            return positions
        positions.append(found)
        start = found + 1


def _same_known_feature(
    existing: dict, motif: CommonMotif, start: int, end: int, direction: int, length: int,
) -> bool:

    existing_name = str(existing.get("name", "")).casefold().replace("×", "x")
    motif_name = motif.name.casefold().replace("×", "x")
    same_name = bool(existing_name) and (motif_name in existing_name or existing_name in motif_name)
    if motif.feature_type == "RBS":
        same_name = same_name or str(existing.get("type", "")).lower() == "rbs" or (
            existing.get("regulatory_class") == "ribosome_binding_site"
        ) or "shine" in existing_name or "sd-like motif" in existing_name
    if existing.get("direction", direction) not in (0, direction):
        return False
    try:
        existing_start = int(existing.get("start", -1))
        existing_end = int(existing.get("end", -1))
    except (TypeError, ValueError):
        return False
    overlap = max(
        max(0, min(existing_end + offset, end) - max(existing_start + offset, start))
        for offset in (-length, 0, length)
    )


    substantially_same_span = overlap / max(1, end - start) >= 0.8
    return same_name and substantially_same_span


def _sd_context(
    dna: str, start: int, end: int, direction: int, circular: bool,
    reverse: str, existing: list[dict],
) -> tuple[str, bool] | None:

    length = len(dna)
    oriented = dna if direction == 1 else reverse
    boundary = end if direction == 1 else length - start
    contexts: list[tuple[str, bool]] = []
    for spacer in range(4, 15):
        codon_start = boundary + spacer
        if circular:
            if end - start + spacer + 3 > length:
                continue
            codon = "".join(oriented[(codon_start + i) % length] for i in range(3))
        else:
            codon = oriented[codon_start:codon_start + 3]
        if codon not in {"ATG", "GTG", "TTG"}:
            continue
        genomic_start = (codon_start if direction == 1 else length - 1 - codon_start) % length
        linked = []
        for feature in existing:
            if str(feature.get("type", "")).lower() != "cds" or feature.get("direction") != direction:
                continue
            first = (
                feature.get("translation_start", feature.get("start")) if direction == 1
                else feature.get("translation_end", feature.get("end", 0)) - 1
            )
            if first is not None and int(first) % length == genomic_start:
                linked.append(str(feature.get("name") or "unnamed CDS"))
        context = (
            f"SD-like motif in the corresponding mRNA, {spacer} nt upstream of a possible "
            f"{codon.replace('T', 'U')} initiation codon (5′→3′); first codon base at DNA "
            f"coordinate {genomic_start + 1} on the {'forward' if direction == 1 else 'reverse'} strand. "
        )
        context += (
            (f"Positionally associated with annotated CDS: {', '.join(name[:200] for name in linked[:3])}"
             + (f" (and {len(linked) - 3} more)" if len(linked) > 3 else '') + '. ') if linked
            else "No annotated CDS starts at this codon; this motif is unassigned. "
        )
        context += "Review required: the 4–14 nt spacer screen does not establish transcription or RBS activity."
        contexts.append((context, bool(linked)))

    return next((context for context in contexts if context[1]), contexts[0] if contexts else None)


def detect_common_features(
    sequence: str,
    *,
    circular: bool = False,
    existing: list[dict] | None = None,
) -> list[dict]:

    dna = clean_dna(sequence)
    if not dna or set(dna) - set("ACGTRYSWKMBDHVN"):
        raise SequenceError("Feature detection requires a non-empty DNA sequence using IUPAC base codes.")
    length = len(dna)
    reverse = reverse_complement(dna)
    existing = existing or []
    matches: list[dict] = []
    seen: set[tuple[str, int, int]] = set()

    for motif in COMMON_MOTIFS:
        motif_dna = validate_dna(motif.sequence, field=f"{motif.name} motif")
        if len(motif_dna) > length:
            continue
        extension = dna[: len(motif_dna) - 1] if circular else ""
        searchable = dna + extension
        strands = ((1, motif_dna), (-1, reverse_complement(motif_dna)))
        for direction, query in strands:
            for start in _occurrences(searchable, query):
                end = start + len(query)
                key = (motif.name, start, end)
                if key in seen:
                    continue
                seen.add(key)
                if any(_same_known_feature(item, motif, start, end, direction, length if circular else 0) for item in existing):
                    continue
                basis = motif.basis
                feature_name, feature_type = motif.name, motif.feature_type
                if motif.feature_type == "RBS":
                    context = _sd_context(dna, start, end, direction, circular, reverse, existing)
                    if context is None:
                        continue
                    basis, linked = context
                    if not linked:
                        feature_name, feature_type = "SD-like motif (unassigned)", "misc_feature"

                    if any(
                        str(item.get("type", "")).lower() == "cds"
                        and item.get("direction") == direction
                        and any(
                            int(item.get("translation_start", item.get("start", 0))) + offset <= start
                            and int(item.get("translation_end", item.get("end", 0))) + offset >= end
                            for offset in (-length, 0, length) if circular or offset == 0
                        ) for item in existing
                    ):
                        continue
                matches.append({
                    "annotation": {
                        "name": feature_name,
                        "type": feature_type,
                        "start": start,
                        "end": end,
                        "direction": direction,
                        "color": motif.color,
                        "inferred": True,
                        "basis": basis,
                        **({"regulatory_class": "ribosome_binding_site"} if feature_type == "RBS" else {}),
                    },
                    "matched_sequence": query,
                    "basis": basis,
                })

    ordered = sorted(
        matches,
        key=lambda match: (
            match["annotation"]["type"] != "RBS",
            match["annotation"]["start"],
            match["annotation"]["end"],
            match["annotation"]["name"],
        ),
    )


    result: list[dict] = []
    for match in ordered:
        feature = match["annotation"]
        if (feature["type"] == "RBS" or feature["name"] == "SD-like motif (unassigned)") and any(
            (previous["annotation"]["type"] == "RBS" or previous["annotation"]["name"] == "SD-like motif (unassigned)")
            and previous["annotation"]["direction"] == feature["direction"]
            and any(
                max(previous["annotation"]["start"] + offset, feature["start"])
                < min(previous["annotation"]["end"] + offset, feature["end"])
                for offset in (-length, 0, length) if circular or offset == 0
            ) for previous in result
        ):
            continue
        result.append(match)
    return sorted(result, key=lambda match: (
        match["annotation"]["start"], match["annotation"]["end"], match["annotation"]["name"],
    ))
