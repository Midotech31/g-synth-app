"""Reproducible design provenance without storing duplicate raw sequences."""
from __future__ import annotations

import hashlib
import json
from typing import Any

from gsynth_engine import __version__
from gsynth_engine.constants import ALL_ENZYMES

PROVENANCE_SCHEMA = "gsynth.provenance/v1"
ENZYME_TABLE_SOURCE = "REBASE geometry via Biopython plus curated G-Synth overrides"


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def enzyme_table_checksum() -> str:
    return sha256_text(_canonical(ALL_ENZYMES))


def _parameter_manifest(parameters: dict[str, Any]) -> dict[str, Any]:
    """Keep choices, but replace large/sensitive molecular inputs by hashes."""
    out: dict[str, Any] = {}
    molecular = {"sequence", "template", "vector_sequence", "insert_reverse", "design"}
    for key, value in parameters.items():
        if key in molecular and isinstance(value, str):
            out[key] = {"length": len(value), "sha256": sha256_text(value.upper())}
        elif key in {"reads", "traces"}:
            out[key] = {"count": len(value) if hasattr(value, "__len__") else None}
        else:
            out[key] = value
    return out


def build_provenance(
    workflow: str,
    *,
    parameters: dict[str, Any],
    output_sequence: str,
    generated_at: str | None = None,
    vector_sequence: str = "",
) -> dict[str, Any]:
    manifest = _parameter_manifest(parameters)
    return {
        "schema": PROVENANCE_SCHEMA,
        "engine_version": __version__,
        "workflow": workflow,
        "generated_at": generated_at,
        "parameters": manifest,
        "parameters_sha256": sha256_text(_canonical(manifest)),
        "output_sha256": sha256_text(output_sequence.upper()),
        "vector_sha256": sha256_text(vector_sequence.upper()) if vector_sequence else None,
        "enzyme_table": {
            "source": ENZYME_TABLE_SOURCE,
            "sha256": enzyme_table_checksum(),
            "entries": len(ALL_ENZYMES),
        },
    }
