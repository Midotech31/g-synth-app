"""SBOL 3 interchange for complete nucleotide constructs."""
from __future__ import annotations

import re

import sbol3

SBOL_NAMESPACE = "https://gsynth.app/designs"

FEATURE_ROLES = {
    "CDS": "https://identifiers.org/SO:0000316",
    "gene": "https://identifiers.org/SO:0000704",
    "promoter": "https://identifiers.org/SO:0000167",
    "terminator": "https://identifiers.org/SO:0000141",
    "RBS": "https://identifiers.org/SO:0000139",
    "rep_origin": "https://identifiers.org/SO:0000296",
}
ROLE_FEATURES = {role: feature_type for feature_type, role in FEATURE_ROLES.items()}


def _display_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]", "_", value.strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_") or "construct"
    return cleaned if cleaned[0].isalpha() else f"construct_{cleaned}"


def to_sbol3(
    sequence: str,
    *,
    name: str = "construct",
    description: str = "",
    features: list[dict] | None = None,
    circular: bool = False,
) -> str:
    """Return a validated SBOL 3 JSON-LD document."""
    display_id = _display_id(name)
    sequence_object = sbol3.Sequence(
        f"{SBOL_NAMESPACE}/{display_id}_sequence",
        elements=sequence.upper(),
        encoding=sbol3.IUPAC_DNA_ENCODING,
        name=f"{name} sequence",
    )
    component = sbol3.Component(
        f"{SBOL_NAMESPACE}/{display_id}",
        [sbol3.SBO_DNA, sbol3.SO_CIRCULAR if circular else sbol3.SO_LINEAR],
        sequences=[sequence_object],
        name=name,
        description=description or None,
    )

    sequence_length = len(sequence)
    for annotation in features or []:
        start = int(annotation.get("start", 0))
        end = int(annotation.get("end", 0))
        if start < 0 or end <= start or start >= sequence_length:
            continue
        orientation = (
            sbol3.SBOL_REVERSE_COMPLEMENT
            if int(annotation.get("direction", 0)) == -1
            else sbol3.SBOL_INLINE
        )
        if end <= sequence_length:
            locations = [sbol3.Range(sequence_object, start + 1, end, orientation=orientation)]
        elif circular:
            locations = [
                sbol3.Range(sequence_object, start + 1, sequence_length, orientation=orientation, order=1),
                sbol3.Range(sequence_object, 1, end - sequence_length, orientation=orientation, order=2),
            ]
        else:
            continue
        role = FEATURE_ROLES.get(str(annotation.get("type", "")))
        component.features.append(sbol3.SequenceFeature(
            locations,
            roles=[role] if role else None,
            orientation=orientation,
            name=str(annotation.get("name") or annotation.get("type") or "feature")[:200],
        ))

    document = sbol3.Document()
    document.add(sequence_object)
    document.add(component)
    report = document.validate()
    if report.errors:
        raise ValueError("Invalid SBOL 3 document: " + "; ".join(report.errors))
    return document.write_string(sbol3.JSONLD)


def read_sbol3(content: str, file_format: str):
    """Read and validate an SBOL 3 document using the reference Python library."""
    formats = {
        "sbol3-jsonld": sbol3.JSONLD,
        "sbol3-turtle": sbol3.TURTLE,
        "sbol3-xml": sbol3.RDF_XML,
    }
    document = sbol3.Document()
    document.read_string(content, formats[file_format])
    report = document.validate()
    if report.errors:
        raise ValueError("; ".join(report.errors))
    return document
