from __future__ import annotations

import textwrap
from dataclasses import dataclass, field

from gsynth_engine.sequence import clean_dna

LINE_WIDTH = 60
GROUP = 10


@dataclass
class Feature:


    name: str
    type: str = "misc_feature"
    start: int = 0
    end: int = 0
    direction: int = 1
    qualifiers: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, entry: dict) -> Feature:
        feature_type = str(entry.get("type") or "misc_feature")
        qualifiers = {}
        regulatory_class = entry.get("regulatory_class") or {
            "RBS": "ribosome_binding_site", "promoter": "promoter", "terminator": "terminator",
        }.get(feature_type)
        if regulatory_class:
            feature_type = "regulatory"
            qualifiers["regulatory_class"] = str(regulatory_class)
        if entry.get("inferred"):
            qualifiers["note"] = "Computational candidate: " + str(entry.get("basis") or "function unconfirmed")
        elif entry.get("basis"):
            qualifiers["note"] = str(entry["basis"])
        start, end = int(entry.get("start", 0)), int(entry.get("end", 0))
        if feature_type == 'CDS' and entry.get('translation_start') is not None and entry.get('translation_end') is not None:


            start, end = int(entry['translation_start']), int(entry['translation_end'])
            qualifiers['codon_start'] = '1'
        return cls(
            name=str(entry.get("name") or entry.get("label") or ""),
            type=feature_type,
            start=start,
            end=end,
            direction=int(entry.get("direction", 1) or 1),
            qualifiers=qualifiers,
        )


def _locus_line(name: str, length: int, *, circular: bool, date: str) -> str:

    identifier = (name or "construct").replace(" ", "_")[:16]
    topology = "circular" if circular else "linear  "
    return (
        f"LOCUS       {identifier:<16} {length:>11} bp    DNA     "
        f"{topology} SYN {date}"
    )


def _location(feature: Feature, length: int) -> str:

    start = feature.start % length if length else feature.start
    end = feature.end

    if length and end > length:

        first = f"{start + 1}..{length}"
        second = f"1..{end - length}"
        span = f"join({first},{second})"
    else:
        span = f"{start + 1}..{end}"

    return f"complement({span})" if feature.direction == -1 else span


def _qualifier(key: str, value: str) -> list[str]:

    escaped = " ".join(str(value).split()).replace('"', '""')
    text = f'/{key}="{escaped}"'
    return textwrap.wrap(
        text, width=79, initial_indent=" " * 21, subsequent_indent=" " * 21,
        break_long_words=False, break_on_hyphens=False,
    ) or [" " * 21 + text]


def to_genbank(
    sequence: str,
    *,
    name: str = "construct",
    description: str = "",
    features: list[dict] | list[Feature] | None = None,
    circular: bool = False,
    date: str = "01-JAN-2000",
    organism: str = "synthetic DNA construct",
    comments: list[str] | None = None,
) -> str:

    seq = clean_dna(sequence).lower()
    length = len(seq)
    entries = [
        item if isinstance(item, Feature) else Feature.from_dict(item)
        for item in (features or [])
    ]

    lines: list[str] = [
        _locus_line(name, length, circular=circular, date=date),
        f"DEFINITION  {description or name}.",
        "ACCESSION   .",
        "VERSION     .",
        "KEYWORDS    .",
        f"SOURCE      {organism}",
        f"  ORGANISM  {organism}",
        "            other sequences; artificial sequences.",
    ]

    for comment in comments or []:
        wrapped = textwrap.wrap(str(comment), width=67) or [""]
        lines.append(f"COMMENT     {wrapped[0]}")
        lines.extend(f"            {line}" for line in wrapped[1:])

    lines.extend([
        "FEATURES             Location/Qualifiers",
        f"     {'source':<16}1..{length}",
        *_qualifier("organism", organism),
        *_qualifier("mol_type", "other DNA"),
    ])

    for feature in entries:
        if feature.end <= feature.start:
            continue
        lines.append(f"     {feature.type[:15]:<16}{_location(feature, length)}")
        if feature.name:
            lines.extend(_qualifier("label", feature.name))
            lines.extend(_qualifier("note", feature.name))
        for key, value in feature.qualifiers.items():
            lines.extend(_qualifier(key, value))

    lines.append("ORIGIN")
    for offset in range(0, length, LINE_WIDTH):
        chunk = seq[offset : offset + LINE_WIDTH]
        groups = " ".join(chunk[i : i + GROUP] for i in range(0, len(chunk), GROUP))
        lines.append(f"{offset + 1:>9} {groups}")
    lines.append("//")

    return "\n".join(lines) + "\n"


def to_fasta(sequence: str, *, name: str = "construct", description: str = "",
             width: int = 70) -> str:

    seq = clean_dna(sequence).upper()
    header = f">{name}" + (f" {description}" if description else "")
    body = [seq[i : i + width] for i in range(0, len(seq), width)] or [""]
    return "\n".join([header, *body]) + "\n"


def oligos_to_fasta(oligos: list[dict], *, key: str = "Name",
                    sequence_key: str = "Sequence (5'->3')") -> str:

    parts = []
    for oligo in oligos:
        parts.append(
            to_fasta(
                str(oligo.get(sequence_key, "")),
                name=str(oligo.get(key, "oligo")),
            )
        )
    return "".join(parts)
