"""Writing GenBank and FASTA — getting a construct out of G-Synth.

A design that only exists inside one program is not finished. Everything
this engine produces — the cassette, the recombinant plasmid, an imported
backbone — comes out as GenBank so it opens in SnapGene, Benchling, ApE or
Geneious with its features intact, and as FASTA when only the bases matter.

**The format is stricter than it looks.** GenBank's LOCUS line is
column-positional, and a parser that cannot read it silently falls back to
defaults: a circular plasmid becomes linear, which changes what a viewer
draws and what a downstream tool computes. The column positions here are
the ones the specification gives, and the tests prove it by writing a
record, parsing it back with an independent parser, and comparing.

**Features keep their strand.** A feature on the bottom strand is written
`complement(start..end)`. Losing that turns a resistance gene into
nonsense on the map, and it is the sort of error that looks fine until
someone tries to use the file.
"""
from __future__ import annotations

import textwrap
from dataclasses import dataclass, field

from gsynth_engine.sequence import clean_dna

#: GenBank wraps sequence at 60 bases, in six groups of ten.
LINE_WIDTH = 60
GROUP = 10


@dataclass
class Feature:
    """One annotation, in the shape the rest of the engine uses."""

    name: str
    type: str = "misc_feature"
    start: int = 0            #: 0-based, half-open
    end: int = 0
    direction: int = 1        #: 1 forward, -1 reverse
    qualifiers: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, entry: dict) -> "Feature":
        return cls(
            name=str(entry.get("name") or entry.get("label") or ""),
            type=str(entry.get("type") or "misc_feature"),
            start=int(entry.get("start", 0)),
            end=int(entry.get("end", 0)),
            direction=int(entry.get("direction", 1) or 1),
        )


def _locus_line(name: str, length: int, *, circular: bool, date: str) -> str:
    """The LOCUS line, in the columns the format actually requires.

    Written by position rather than by joining with spaces: a parser reads
    `circular` from columns 56–63, and a line that merely looks right gives
    a linear record with no error message.
    """
    identifier = (name or "construct").replace(" ", "_")[:16]
    topology = "circular" if circular else "linear  "
    return (
        f"LOCUS       {identifier:<16} {length:>11} bp    DNA     "
        f"{topology} SYN {date}"
    )


def _location(feature: Feature, length: int) -> str:
    """A GenBank location, wrapping and strand included."""
    start = feature.start % length if length else feature.start
    end = feature.end

    if length and end > length:
        # Straddles the origin: GenBank spells that as a join.
        first = f"{start + 1}..{length}"
        second = f"1..{end - length}"
        span = f"join({first},{second})"
    else:
        span = f"{start + 1}..{end}"

    return f"complement({span})" if feature.direction == -1 else span


def _qualifier(key: str, value: str) -> list[str]:
    """A wrapped /key="value" block at GenBank's 21-column indent."""
    text = f'/{key}="{value}"'
    return textwrap.wrap(
        text, width=58, initial_indent=" " * 21, subsequent_indent=" " * 21
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
) -> str:
    """Render a sequence and its features as a GenBank record.

    Args:
        date: supplied rather than taken from the clock, so the same design
            always produces the same file. Callers that want today's date
            pass it in.

    Returns:
        A complete record, terminated with `//`.
    """
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
        "FEATURES             Location/Qualifiers",
        f"     {'source':<16}1..{length}",
        *_qualifier("organism", organism),
        *_qualifier("mol_type", "other DNA"),
    ]

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
    """Render a sequence as FASTA. Bases only — features do not survive."""
    seq = clean_dna(sequence).upper()
    header = f">{name}" + (f" {description}" if description else "")
    body = [seq[i : i + width] for i in range(0, len(seq), width)] or [""]
    return "\n".join([header, *body]) + "\n"


def oligos_to_fasta(oligos: list[dict], *, key: str = "Name",
                    sequence_key: str = "Sequence (5'->3')") -> str:
    """Every oligo in an order sheet as one FASTA file.

    Suppliers accept a FASTA upload; retyping thirty oligo names into a web
    form is where transcription errors come from.
    """
    parts = []
    for oligo in oligos:
        parts.append(
            to_fasta(
                str(oligo.get(sequence_key, "")),
                name=str(oligo.get(key, "oligo")),
            )
        )
    return "".join(parts)
