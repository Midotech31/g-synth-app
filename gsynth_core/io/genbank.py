"""
GenBank format reader and writer — pure Python, no Biopython dependency.

Handles:
- Multi-record files.
- LOCUS / DEFINITION / ACCESSION / VERSION / KEYWORDS headers.
- FEATURES table (CDS, gene, regulatory, etc.) with qualifiers.
- Compound locations (join, complement, order).
- ORIGIN sequence block.
- Round-trips for the subset of features G-Synth uses.

For very complex or non-standard GenBank files, fall back to Biopython
(`Bio.SeqIO.parse(..., "genbank")`).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, TextIO

from gsynth_core.io.records import Feature, FeatureLocation, SeqRecord
from gsynth_core.sequence.ops import clean_dna


_LOCUS_RE = re.compile(
    r"LOCUS\s+(?P<name>\S+)\s+(?P<length>\d+)\s+bp\s+(?P<mol>\S+)\s+"
    r"(?P<topo>linear|circular)?",
    re.IGNORECASE,
)
_FEATURE_HEADER_RE = re.compile(r"^\s{5}(\S+)\s+(.+)$")
_QUALIFIER_RE = re.compile(r"^\s{21}/(\S+?)=(.*)$|^\s{21}/(\S+)$")
_LOCATION_RANGE_RE = re.compile(r"(\d+)\.\.(\d+)")
_LOCATION_SINGLE_RE = re.compile(r"^(\d+)$")


def _parse_location(loc_str: str) -> FeatureLocation:
    """Parse a GenBank location string into a :class:`FeatureLocation`."""
    s = loc_str.strip()
    strand = 1
    if s.startswith("complement("):
        strand = -1
        s = s[len("complement(") : -1]
    if s.startswith("join(") or s.startswith("order("):
        prefix_len = 5
        inner = s[prefix_len : -1]
        parts: list[tuple[int, int]] = []
        for piece in inner.split(","):
            piece = piece.strip()
            m = _LOCATION_RANGE_RE.fullmatch(piece)
            if m:
                parts.append((int(m.group(1)) - 1, int(m.group(2))))
            else:
                m2 = _LOCATION_SINGLE_RE.match(piece)
                if m2:
                    p = int(m2.group(1)) - 1
                    parts.append((p, p + 1))
        if not parts:
            return FeatureLocation(start=0, end=0, strand=strand)
        return FeatureLocation(
            start=min(p[0] for p in parts),
            end=max(p[1] for p in parts),
            strand=strand,
            compound_parts=tuple(parts),
        )
    m = _LOCATION_RANGE_RE.fullmatch(s)
    if m:
        return FeatureLocation(start=int(m.group(1)) - 1, end=int(m.group(2)), strand=strand)
    m2 = _LOCATION_SINGLE_RE.match(s)
    if m2:
        p = int(m2.group(1)) - 1
        return FeatureLocation(start=p, end=p + 1, strand=strand)
    return FeatureLocation(start=0, end=0, strand=strand)


def _parse_genbank_stream(stream: Iterable[str]) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    name = ""
    description = ""
    topology: str = "linear"
    molecule_type: str = "DNA"
    annotations: dict[str, str] = {}
    features: list[Feature] = []
    sequence_lines: list[str] = []

    state: str = "header"
    current_feature_type = ""
    current_location_parts: list[str] = []
    current_qualifiers: dict[str, list[str]] = {}
    current_qual_key = ""

    def finalize_current_feature() -> None:
        nonlocal current_feature_type, current_location_parts, current_qualifiers
        if not current_feature_type:
            return
        loc = _parse_location("".join(current_location_parts))
        qmap: dict[str, str] = {}
        for k, v in current_qualifiers.items():
            qmap[k] = "\n".join(v).strip().strip('"')
        features.append(Feature(type=current_feature_type, location=loc, qualifiers=qmap))
        current_feature_type = ""
        current_location_parts = []
        current_qualifiers = {}

    def finalize_current_record() -> None:
        nonlocal name, description, topology, molecule_type, annotations
        nonlocal features, sequence_lines, state
        finalize_current_feature()
        seq = clean_dna("".join(sequence_lines), allow_ambiguous=True)
        if name or seq:
            records.append(SeqRecord(
                sequence=seq,
                name=name,
                description=description,
                features=features,
                topology="circular" if topology == "circular" else "linear",
                molecule_type=molecule_type,
                annotations=dict(annotations),
            ))
        name = ""
        description = ""
        topology = "linear"
        molecule_type = "DNA"
        annotations = {}
        features = []
        sequence_lines = []
        state = "header"

    for raw in stream:
        line = raw.rstrip("\r\n")

        if line.startswith("LOCUS"):
            finalize_current_record()
            m = _LOCUS_RE.match(line)
            if m:
                name = m.group("name")
                molecule_type = m.group("mol")
                topology = (m.group("topo") or "linear").lower()
            continue

        if line.startswith("DEFINITION"):
            description = line[12:].strip()
            state = "header"
            continue

        if line.startswith("ACCESSION"):
            annotations["accession"] = line[12:].strip()
            continue
        if line.startswith("VERSION"):
            annotations["version"] = line[12:].strip()
            continue
        if line.startswith("KEYWORDS"):
            annotations["keywords"] = line[12:].strip()
            continue
        if line.startswith("SOURCE"):
            annotations["source"] = line[12:].strip()
            continue

        if line.startswith("FEATURES"):
            state = "features"
            continue

        if line.startswith("ORIGIN"):
            finalize_current_feature()
            state = "origin"
            continue

        if line.startswith("//"):
            finalize_current_record()
            continue

        if state == "features":
            # Feature header: 5 leading spaces, then type, then location
            m = _FEATURE_HEADER_RE.match(line)
            if m:
                finalize_current_feature()
                current_feature_type = m.group(1)
                current_location_parts = [m.group(2).strip()]
                current_qualifiers = {}
                current_qual_key = ""
                continue
            # Qualifier line: 21 leading spaces, then /key=value
            mq = _QUALIFIER_RE.match(line)
            if mq:
                if mq.group(1):
                    current_qual_key = mq.group(1)
                    current_qualifiers.setdefault(current_qual_key, []).append(mq.group(2))
                else:
                    current_qual_key = mq.group(3)
                    current_qualifiers.setdefault(current_qual_key, []).append("")
                continue
            # Continuation of location (no slash on this line)
            stripped = line.strip()
            if stripped and current_feature_type and not stripped.startswith("/"):
                if current_qual_key:
                    current_qualifiers[current_qual_key][-1] += " " + stripped
                else:
                    current_location_parts.append(stripped)
                continue

        elif state == "origin":
            # Origin lines: "        1 atgaaacgcc tggctgttt t ..."
            digits_stripped = re.sub(r"[\s\d]", "", line)
            if digits_stripped:
                sequence_lines.append(digits_stripped.upper())

    finalize_current_record()
    return records


def parse_genbank(source: str | Path | TextIO) -> list[SeqRecord]:
    """Parse one or more GenBank records from file or string."""
    if isinstance(source, (str, Path)):
        p = Path(source)
        if p.exists() and p.is_file():
            with p.open("r", encoding="utf-8") as f:
                return _parse_genbank_stream(f)
        if isinstance(source, str) and "LOCUS" in source:
            return _parse_genbank_stream(source.splitlines(keepends=True))
        raise FileNotFoundError(f"Not a valid GenBank file or string: {source!r}")
    return _parse_genbank_stream(source)


def _format_location(loc: FeatureLocation) -> str:
    if loc.compound_parts:
        body = "join(" + ",".join(f"{s+1}..{e}" for s, e in loc.compound_parts) + ")"
    else:
        body = f"{loc.start + 1}..{loc.end}"
    if loc.strand == -1:
        body = f"complement({body})"
    return body


def _format_qualifier(key: str, value: str) -> list[str]:
    """Wrap a qualifier line at 80 chars (column 21 indent)."""
    indent = " " * 21
    if not value:
        body = f"/{key}"
    elif value.startswith('"') or value.lstrip("-").isdigit():
        body = f"/{key}={value}"
    else:
        body = f'/{key}="{value}"'
    width = 58  # 80 − 21 indent − 1 for slash
    if len(body) <= width:
        return [indent + body]
    out: list[str] = []
    while body:
        out.append(indent + body[:width])
        body = body[width:]
    return out


def write_genbank(
    records: SeqRecord | list[SeqRecord],
    destination: str | Path | TextIO | None = None,
) -> str:
    """Serialize records to GenBank text. Returns the string."""
    if isinstance(records, SeqRecord):
        records = [records]
    chunks: list[str] = []
    for rec in records:
        L = len(rec.sequence)
        topo = "circular" if rec.topology == "circular" else "linear  "
        chunks.append(
            f"LOCUS       {rec.name or 'unnamed':<16s} {L:>10d} bp    "
            f"{rec.molecule_type:<7s} {topo:<8s}"
        )
        if rec.description:
            chunks.append(f"DEFINITION  {rec.description}")
        for key in ("accession", "version", "keywords", "source"):
            if key in rec.annotations:
                chunks.append(f"{key.upper():<12s}{rec.annotations[key]}")
        chunks.append("FEATURES             Location/Qualifiers")
        for f in rec.features:
            chunks.append(f"     {f.type:<16s}{_format_location(f.location)}")
            for k, v in f.qualifiers.items():
                chunks.extend(_format_qualifier(k, str(v)))
        chunks.append("ORIGIN")
        seq = rec.sequence.lower()
        for i in range(0, len(seq), 60):
            line = seq[i : i + 60]
            blocks = " ".join(line[j : j + 10] for j in range(0, len(line), 10))
            chunks.append(f"{i + 1:>9d} {blocks}")
        chunks.append("//")

    out = "\n".join(chunks) + "\n"

    if destination is not None:
        if isinstance(destination, (str, Path)):
            Path(destination).write_text(out, encoding="utf-8")
        else:
            destination.write(out)
    return out
