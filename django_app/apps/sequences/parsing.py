from __future__ import annotations

import io
import re
from dataclasses import asdict, dataclass, field

import sbol3
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

from apps.sequences.sbol import ROLE_FEATURES, feature_metadata, read_sbol3

MAX_UPLOAD_BYTES = 10 * 1024 * 1024


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


SKIPPED_FEATURE_TYPES = frozenset({"source"})


class ParseError(ValueError):
    pass


@dataclass
class Annotation:

    name: str
    type: str
    start: int
    end: int
    direction: int
    color: str
    regulatory_class: str = ""
    inferred: bool = False
    basis: str = ""
    translation_start: int | None = None
    translation_end: int | None = None


@dataclass
class ParsedRecord:
    name: str
    description: str
    sequence: str
    length: int
    topology: str
    gc_content: float
    source_format: str
    annotations: list[Annotation] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {**asdict(self),
                "annotations": [{key: value for key, value in asdict(a).items()
                                 if value is not None and (key not in {"regulatory_class", "basis", "inferred"} or value)}
                                for a in self.annotations]}


SNAPGENE_MAGIC = b"\x09\x00\x00\x00\x0e"


def detect_format(text: str | bytes, filename: str = "") -> str:

    lower = filename.lower()

    if isinstance(text, bytes):
        if text.startswith(SNAPGENE_MAGIC) or lower.endswith((".dna", ".prot")):
            return "snapgene"
        try:
            text = text.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise ParseError(
                "This file is not text and is not a SnapGene file. Expected "
                "FASTA, GenBank, or SnapGene's own .dna format."
            ) from error

    head = text.lstrip()[:200].upper()
    sbol_marker = "http://sbols.org/v3#"
    if lower.endswith((".jsonld", ".sbol.json")) or (
        text.lstrip().startswith("{") and sbol_marker in text
    ):
        return "sbol3-jsonld"
    if lower.endswith((".ttl", ".sbol")) or (
        sbol_marker in text and ("@prefix" in text.lower() or "PREFIX" in text[:200])
    ):
        return "sbol3-turtle"
    if lower.endswith((".rdf", ".xml")) and sbol_marker in text:
        return "sbol3-xml"
    if head.startswith("LOCUS"):
        return "genbank"
    if head.startswith(">"):
        return "fasta"
    if lower.endswith((".gb", ".gbk", ".genbank", ".ape")):
        return "genbank"
    if lower.endswith((".fa", ".fasta", ".fna", ".ffn", ".faa", ".seq")):
        return "fasta"
    raise ParseError(
        "Unrecognised file. Expected FASTA (starting with '>'), GenBank "
        "(starting with 'LOCUS'), SBOL 3, or a SnapGene .dna file."
    )


def _record_from_sbol(content: str, source_format: str, filename: str) -> ParsedRecord:
    document = read_sbol3(content, source_format)
    components = [item for item in document if isinstance(item, sbol3.Component)]
    component = next((item for item in components if item.sequences), None)
    if component is None:
        raise ParseError("The SBOL 3 document contains no DNA component with a sequence.")
    sequence_object = document.find(component.sequences[0])
    if not isinstance(sequence_object, sbol3.Sequence) or not sequence_object.elements:
        raise ParseError("The SBOL 3 component does not contain nucleotide elements.")

    sequence = sequence_object.elements.upper()
    annotations: list[Annotation] = []
    for feature in component.features:
        if not isinstance(feature, sbol3.SequenceFeature):
            continue
        ranges = [location for location in feature.locations if isinstance(location, sbol3.Range)]
        if not ranges:
            continue
        if any(str(location.sequence) != sequence_object.identity or not 1 <= location.start <= location.end <= len(sequence) for location in ranges):
            raise ParseError('SBOL feature locations must refer to the selected sequence within its bounds.')
        orientations = {location.orientation or feature.orientation for location in ranges}
        if len(orientations) != 1:
            raise ParseError('Mixed-strand SBOL features cannot be displayed as one contiguous span.')
        start, end = ranges[0].start - 1, ranges[0].end
        if len(ranges) > 1:
            boundary = next((part for part in ranges if part.end == len(sequence)), None)
            origin = next((part for part in ranges if part.start == 1), None)
            if (sbol3.SO_CIRCULAR not in component.types or len(ranges) != 2 or boundary is None or origin is None
                    or boundary is origin or origin.end >= boundary.start):
                raise ParseError('Discontinuous SBOL features cannot be displayed as one contiguous span.')
            start, end = boundary.start - 1, len(sequence) + origin.end
        metadata = feature_metadata(feature)
        role = next((value for value in feature.roles if value in ROLE_FEATURES), "")
        feature_type = str(metadata.get('type') or ROLE_FEATURES.get(role, "misc_feature"))
        orientation = ranges[0].orientation or feature.orientation
        color = str(metadata.get('color') or FEATURE_COLORS.get(feature_type, DEFAULT_FEATURE_COLOR))
        if not re.fullmatch(r'#[0-9a-fA-F]{6}', color):
            raise ParseError('SBOL feature color must be a six-digit hexadecimal color.')
        translation_start, translation_end = metadata.get('translation_start'), metadata.get('translation_end')
        if translation_start is not None or translation_end is not None:
            if (type(translation_start) is not int or type(translation_end) is not int
                    or not start <= translation_start < translation_end <= end):
                raise ParseError('SBOL coding bounds must lie inside the feature.')
        if type(metadata.get('inferred', False)) is not bool:
            raise ParseError('SBOL candidate status must be boolean.')
        annotations.append(Annotation(
            name=(feature.name or feature.display_id or feature_type)[:200],
            type=feature_type,
            start=start,
            end=end,
            direction=-1 if orientation == sbol3.SBOL_REVERSE_COMPLEMENT else 1 if orientation == sbol3.SBOL_INLINE else 0,
            color=color,
            inferred=metadata.get('inferred', False),
            regulatory_class=str(metadata.get('regulatory_class') or ''),
            basis=(feature.description or '')[:4000],
            translation_start=translation_start,
            translation_end=translation_end,
        ))
    annotations.sort(key=lambda annotation: (annotation.start, -(annotation.end - annotation.start)))

    topology = "circular" if sbol3.SO_CIRCULAR in component.types else "linear"
    fallback = filename.rsplit("/", 1)[-1].split(".", 1)[0] or "Untitled"
    return ParsedRecord(
        name=(component.name or component.display_id or fallback)[:120],
        description=(component.description or "")[:300],
        sequence=sequence,
        length=len(sequence),
        topology=topology,
        gc_content=_gc_content(sequence),
        source_format="sbol3",
        annotations=annotations,
    )


def _feature_name(feature: SeqFeature, fallback: str) -> str:

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
            continue


        parts = list(getattr(location, "parts", ()))
        record_length = len(record.seq)
        boundary_part = next(
            (part for part in parts if int(part.end) == record_length), None
        )
        origin_part = next((part for part in parts if int(part.start) == 0), None)
        is_circular = str(record.annotations.get("topology", "")).lower() == "circular"
        if (
            is_circular
            and boundary_part is not None
            and origin_part is not None
            and len(parts) == 2
        ):
            start = int(boundary_part.start)
            end = record_length + int(origin_part.end)
        if len(parts) > 1 and (
            not is_circular or len(parts) != 2 or boundary_part is None or origin_part is None
            or sum(len(part) for part in parts) != end - start
        ):
            raise ParseError(
                f"Feature '{_feature_name(feature, ftype)}' has a discontinuous location. "
                "G-Synth cannot represent this feature faithfully as one contiguous span."
            )
        if end <= start:
            continue
        translation_start = translation_end = None
        if ftype == 'CDS':
            if str(feature.qualifiers.get('transl_table', ['1'])[0]) not in {'1', '11'} or feature.qualifiers.get('transl_except'):
                raise ParseError('This CDS uses an unsupported genetic code or translation exception; its protein cannot be displayed faithfully.')
            try:
                offset = int(feature.qualifiers.get('codon_start', ['1'])[0]) - 1
            except (ValueError, TypeError):
                raise ParseError('CDS codon_start must be 1, 2 or 3.') from None
            if offset not in (0, 1, 2) or end - start <= offset:
                raise ParseError('CDS codon_start must identify a base inside the feature.')
            translation_start = start + (offset if location.strand != -1 else 0)
            translation_end = end - (offset if location.strand == -1 else 0)
        out.append(Annotation(
            name=_feature_name(feature, ftype),
            type=ftype,
            start=start,
            end=end,
            direction=(
                1 if location.strand == 1 else -1 if location.strand == -1 else 0
            ),
            color=FEATURE_COLORS.get(ftype, DEFAULT_FEATURE_COLOR),
            translation_start=translation_start,
            translation_end=translation_end,
            regulatory_class=str(feature.qualifiers.get("regulatory_class", [""])[0]),
            inferred=bool(feature.qualifiers.get("inference")) or any(
                note.startswith("Computational candidate:") for note in feature.qualifiers.get("note", [])
            ),
            basis=next((note.removeprefix("Computational candidate: ") for note in feature.qualifiers.get("note", [])
                        if note.startswith("Computational candidate:")),
                       ' '.join(str(note) for note in feature.qualifiers.get('note', [])
                                if str(note) != _feature_name(feature, ftype))
                       or ' '.join(str(value) for value in feature.qualifiers.get('inference', [])))[:4000],
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

    return "linear"


def parse_sequence_file(content: bytes | str, filename: str = "") -> ParsedRecord:

    if isinstance(content, bytes) and len(content) > MAX_UPLOAD_BYTES:
        raise ParseError(
            f"File is larger than {MAX_UPLOAD_BYTES // (1024 * 1024)} MB."
        )


    if not content or not (
        content.strip() if isinstance(content, str) else content.strip()
    ):
        raise ParseError("The file is empty.")

    source_format = detect_format(content, filename)

    if source_format.startswith("sbol3-"):
        if isinstance(content, bytes):
            content = content.decode("utf-8")
        try:
            return _record_from_sbol(content, source_format, filename)
        except Exception as exc:
            if isinstance(exc, ParseError):
                raise
            raise ParseError(f"Could not read this SBOL 3 file: {exc}") from exc

    if source_format == "snapgene":
        if isinstance(content, str):
            raise ParseError("A SnapGene file must be uploaded as a file, not text.")
        handle: io.IOBase = io.BytesIO(content)
    else:
        if isinstance(content, bytes):
            try:
                text = content.decode("utf-8")
            except UnicodeDecodeError:
                text = content.decode("latin-1", errors="replace")
        else:
            text = content
        if not text.strip():
            raise ParseError("The file is empty.")
        handle = io.StringIO(text)

    try:
        records = list(SeqIO.parse(handle, source_format))
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

        annotations=_annotations_from(record),
    )
