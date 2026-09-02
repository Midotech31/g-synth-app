"""Generate machine-readable evidence for peptide starts and enzyme coverage."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from gsynth_engine.cloning import translate
from gsynth_engine.codon import TABLES, optimise
from gsynth_engine.constants import ALL_ENZYMES, overhang

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "publication_evidence" / "peptide_and_enzyme_validation.json"
MATURE_GLARGINE_A = "GIVEQCCTSICSLYQLENYCG"


def digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def peptide_case(sequence: str, context: str) -> dict:
    result = optimise(
        sequence,
        is_protein=True,
        protein_context=context,
        keep_stop=False,
    )
    translated = translate(result.sequence)
    return {
        "input_protein": result.input_protein,
        "requested_context": context,
        "resolved_context": result.protein_context,
        "encoded_protein": result.protein,
        "initiator_methionine_added": result.initiator_methionine_added,
        "recommended_design_is_coding": result.recommended_design_is_coding,
        "translation_matches_encoded_protein": translated == result.protein,
        "dna_length_nt": result.length,
        "dna_sha256": digest(result.sequence),
    }


def main() -> None:
    aliases = [
        alias
        for spec in ALL_ENZYMES.values()
        for alias in spec.get("aliases", ())
    ]
    host_results = {}
    for key, table in TABLES.items():
        result = optimise(
            MATURE_GLARGINE_A,
            table=table,
            is_protein=True,
            protein_context="auto",
            keep_stop=False,
        )
        host_results[key] = {
            "translated_protein": translate(result.sequence),
            "protein_preserved": translate(result.sequence) == MATURE_GLARGINE_A,
            "starts_with_atg": result.sequence.startswith("ATG"),
            "dna_sha256": digest(result.sequence),
        }

    hindiii_overhang, hindiii_polarity = overhang("HindIII")
    evidence = {
        "schema_version": "1.0",
        "peptide_cases": {
            "automatic_mature_glargine_a": peptide_case(MATURE_GLARGINE_A, "auto"),
            "automatic_existing_methionine": peptide_case("MGIVEQ", "auto"),
            "explicit_complete_orf_without_methionine": peptide_case(
                "GIVEQ", "complete_orf",
            ),
            "explicit_mature_peptide_starting_with_methionine": peptide_case(
                "MGIVEQ", "mature_peptide",
            ),
        },
        "host_back_translation": {
            "host_count": len(host_results),
            "all_preserve_mature_glargine_a": all(
                row["protein_preserved"] for row in host_results.values()
            ),
            "all_avoid_artificial_atg": all(
                not row["starts_with_atg"] for row in host_results.values()
            ),
            "hosts": host_results,
        },
        "restriction_catalogue": {
            "canonical_cut_geometries": len(ALL_ENZYMES),
            "selectable_names": len(ALL_ENZYMES) + len(aliases),
            "alias_names_are_unique": len(aliases) == len(set(aliases)),
            "hindiii": {
                "recognition": ALL_ENZYMES["HindIII"]["recognition"],
                "overhang": hindiii_overhang,
                "polarity": hindiii_polarity,
            },
        },
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
