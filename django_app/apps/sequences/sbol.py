from __future__ import annotations

import json
import re

import sbol3

SBOL_NAMESPACE = "https://gsynth.app/designs"
ANNOTATION_PROPERTY = "https://gsynth.app/terms/annotation"
ANNOTATION_FIELDS = {"type", "color", "regulatory_class", "inferred", "translation_start", "translation_end"}


def feature_metadata(feature):
    raw = sbol3.TextProperty(feature, ANNOTATION_PROPERTY, 0, 1).get()
    metadata = json.loads(raw) if raw else {}
    if not isinstance(metadata, dict):
        raise ValueError("SBOL annotation metadata must be an object.")
    return {key: value for key, value in metadata.items() if key in ANNOTATION_FIELDS}

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
        if not 0 <= start < sequence_length or not start < end <= start + sequence_length:
            raise ValueError("SBOL feature coordinates are outside the sequence.")
        direction = int(annotation.get("direction", 0))
        orientation = (
            sbol3.SBOL_REVERSE_COMPLEMENT
            if direction == -1
            else sbol3.SBOL_INLINE if direction == 1 else None
        )
        if end <= sequence_length:
            locations = [sbol3.Range(sequence_object, start + 1, end, orientation=orientation)]
        elif circular:
            locations = [
                sbol3.Range(sequence_object, start + 1, sequence_length, orientation=orientation, order=2 if direction == -1 else 1),
                sbol3.Range(sequence_object, 1, end - sequence_length, orientation=orientation, order=1 if direction == -1 else 2),
            ]
        else:
            raise ValueError("A linear SBOL feature cannot cross the sequence boundary.")
        role = FEATURE_ROLES.get(str(annotation.get("type", "")))
        feature = sbol3.SequenceFeature(
            locations,
            roles=[role] if role else None,
            orientation=orientation,
            name=str(annotation.get("name") or annotation.get("type") or "feature")[:200],
            description=str(annotation.get("basis") or "")[:4000] or None,
        )
        feature.gsynth_annotation = sbol3.TextProperty(feature, ANNOTATION_PROPERTY, 0, 1,
            initial_value=json.dumps({key: annotation[key] for key in ANNOTATION_FIELDS if key in annotation}, sort_keys=True))
        component.features.append(feature)

    document = sbol3.Document()
    document.add(sequence_object)
    document.add(component)
    report = document.validate()
    if report.errors:
        raise ValueError("Invalid SBOL 3 document: " + "; ".join(report.errors))
    return document.write_string(sbol3.JSONLD)


def read_sbol3(content: str, file_format: str):

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
