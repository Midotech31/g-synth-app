"""Sequence file parsing — FASTA and GenBank in, viewer-ready JSON out.

The browser receives annotations already shaped for SeqViz
(`name/start/end/direction/color`), so the frontend does no biology and no
coordinate arithmetic — the two places this kind of code usually goes wrong.

Coordinates are 0-based, half-open (Python/Biopython convention), which is
also what SeqViz expects.
"""
from __future__ import annotations

import io
from dataclasses import asdict, dataclass, field

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

MAX_UPLOAD_BYTES = 10 * 1024 * 1024   # 10 MB

# Feature types we surface, with the colours the viewer draws them in.
# Anything not listed is still shown, in the neutral fallback colour.
FEATURE_COLORS: dict[str, str] = {
    "CDS":           "#0E6E77",
    "gene":          "#2A9D8F",
    "promoter":      "#C97634",
    "terminator":    "#9E3D3D",
    "rep_origin":    "#6A4C93",
    "oriT":          "#6A4C93",
    "regulatory":    "#B8860B",
    "RBS":           "#B8860B",
    "misc_feature":  "#5D6B7A",
    "primer_bind":   "#4A7C59",
    "protein_bind":  "#4A7C59",
    "source":        "#A9B4C0",
    "tRNA":          "#D08C60",
    "rRNA":          "#D08C60",
}
DEFAULT_FEATURE_COLOR = "#5D6B7A"

# Features that carry no useful visual information for a plasmid map.
SKIPPED_FEATURE_TYPES = frozenset({"source"})


class ParseError(ValueError):
    """Raised when a file cannot be read as FASTA or GenBank."""


@dataclass
class Annotation:
    """One drawable feature, already in SeqViz's shape."""
    name: str
    type: str
    start: int
    end: int
    direction: int          # 1 forward, -1 reverse, 0 unstranded
    color: str


@dataclass
class ParsedRecord:
    name: str
    description: str
    sequence: str
    length: int
    topology: str           # "circular" | "linear"
    gc_content: float       # percent, 1 decimal
    source_format: str      # "fasta" | "genbank"
    annotations: list[Annotation] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {**asdict(self),
                "annotations": [asdict(a) for a in self.annotations]}


def detect_format(text: str, filename: str = "") -> str:
    """Sniff the format from the content first, filename second."""
    head = text.lstrip()[:200].upper()
    if head.startswith("LOCUS"):
        return "genbank"
    if head.startswith(">"):
        return "fasta"
    lower = filename.lower()
    if lower.endswith((".gb", ".gbk", ".genbank", ".ape")):
        return "genbank"
    if lower.endswith((".fa", ".fasta", ".fna", ".ffn", ".faa", ".seq")):
        return "fasta"
    raise ParseError(
        "Unrecognised file. Expected FASTA (starting with '>') or "
        "GenBank (starting with 'LOCUS')."
    )


def _feature_name(feature: SeqFeature, fallback: str) -> str:
    """Best human label for a feature, in the order a biologist would pick."""
    for key in ("label", "gene", "product", "note", "standard_name",
                "bound_moiety", "organism"):
        values = feature.qualifiers.get(key)
        if values:
            name = str(values[0]).strip()
            if name:
                return name[:80]
    return fallback


def _annotations_from(record) -> list[Annotation]:
    out: list[Annotation] = []
    for feature in record.features:
        ftype = str(feature.type)
        if ftype in SKIPPED_FEATURE_TYPES:
            continue
        location = feature.location
        if location is None:
            continue
        try:
            start, end = int(location.start), int(location.end)
        except (TypeError, ValueError):
            continue      # fuzzy/unknown locations — nothing sensible to draw
        if end <= start:
            continue
        out.append(Annotation(
            name=_feature_name(feature, ftype),
            type=ftype,
            start=start,
            end=end,
            direction=1 if location.strand in (None, 1) else -1,
            color=FEATURE_COLORS.get(ftype, DEFAULT_FEATURE_COLOR),
        ))
    out.sort(key=lambda a: (a.start, -(a.end - a.start)))
    return out


def _gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    counted = sum(1 for base in sequence.upper() if base in "ACGT")
    if not counted:
        return 0.0
    gc = sum(1 for base in sequence.upper() if base in "GC")
    return round(100.0 * gc / counted, 1)


def _topology(record, source_format: str) -> str:
    topology = str(record.annotations.get("topology", "")).lower()
    if topology in ("circular", "linear"):
        return topology
    # FASTA carries no topology; assume linear rather than guessing.
    return "linear"


def parse_sequence_file(content: bytes | str, filename: str = "") -> ParsedRecord:
    """Parse the first record of a FASTA or GenBank file.

    Raises ParseError with a message meant for an end user, not a stack trace.
    """
    if isinstance(content, bytes):
        if len(content) > MAX_UPLOAD_BYTES:
            raise ParseError(
                f"File is larger than {MAX_UPLOAD_BYTES // (1024 * 1024)} MB."
            )
        try:
            text = content.decode("utf-8")
        except UnicodeDecodeError:
            text = content.decode("latin-1", errors="replace")
    else:
        text = content

    if not text.strip():
        raise ParseError("The file is empty.")

    source_format = detect_format(text, filename)

    try:
        records = list(SeqIO.parse(io.StringIO(text), source_format))
    except Exception as exc:                       # noqa: BLE001
        raise ParseError(f"Could not read this {source_format} file: {exc}") from exc

    if not records:
        raise ParseError(f"No sequence records found in this {source_format} file.")

    record = records[0]
    sequence = str(record.seq).upper()
    if not sequence:
        raise ParseError("The record contains no sequence.")

    name = (record.name or record.id or "").strip()
    if name in ("", "<unknown name>", "<unknown id>"):
        name = (filename.rsplit("/", 1)[-1].rsplit(".", 1)[0] or "Untitled")

    description = (record.description or "").strip()
    if description == name:
        description = ""

    return ParsedRecord(
        name=name[:120],
        description=description[:300],
        sequence=sequence,
        length=len(sequence),
        topology=_topology(record, source_format),
        gc_content=_gc_content(sequence),
        source_format=source_format,
        annotations=_annotations_from(record) if source_format == "genbank" else [],
    )
