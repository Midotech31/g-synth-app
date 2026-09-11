#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

from gsynth_engine.cloning import translate
from gsynth_engine.codon import (
    CODON_DATA_SHA256,
    CODON_DATA_VERSION,
    SYNONYMS,
    TABLES,
    back_translate,
)

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "publication_evidence" / "codon_host_profile_validation.json"
BENCHMARK_PROTEIN = "ACDEFGHIKLMNPQRSTVWY" * 5


def _weight_signature(weights: dict[str, float]) -> str:
    payload = "\n".join(
        f"{codon}\t{weights[codon]:.15g}" for codon in sorted(weights)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def main() -> None:
    records = []
    generated_sequences: set[str] = set()
    weight_signatures: set[str] = set()
    for key, table in TABLES.items():
        sequence = back_translate(BENCHMARK_PROTEIN, table=table)
        signature = _weight_signature(table.weights)
        generated_sequences.add(sequence)
        weight_signatures.add(signature)
        records.append({
            "key": key,
            "name": table.name,
            "category": table.category,
            "ncbi_taxon_id": table.taxon_id,
            "dataset": table.dataset,
            "dataset_release": table.dataset_release,
            "data_scope": table.data_scope,
            "coding_sequences": table.coding_sequences,
            "codon_count": table.codon_count,
            "gc_percent": table.gc_percent,
            "codons_present": len(table.weights),
            "synonymous_families_normalized": all(
                max(table.weight(codon) for codon in codons) == 1.0
                for codons in SYNONYMS.values()
            ),
            "weight_vector_sha256": signature,
            "benchmark_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "benchmark_translation_preserved": translate(sequence) == BENCHMARK_PROTEIN,
        })

    evidence = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "source_snapshot": {
            "dataset": "FDA HIVE-CUTs / CoCoPUTs",
            "release": CODON_DATA_VERSION,
            "sha256": CODON_DATA_SHA256,
            "method": (
                "Raw codon counts normalized to the maximum count within each "
                "synonymous amino-acid family"
            ),
        },
        "validation": {
            "host_profiles": len(TABLES),
            "complete_64_codon_profiles": sum(
                len(table.weights) == 64 for table in TABLES.values()
            ),
            "distinct_weight_vectors": len(weight_signatures),
            "distinct_preferred_codon_benchmark_sequences": len(generated_sequences),
            "benchmark_length_residues": len(BENCHMARK_PROTEIN),
            "protein_invariant_passed_for_all_hosts": all(
                record["benchmark_translation_preserved"] for record in records
            ),
        },
        "interpretation": {
            "metric": "Profile-relative CAI for bundled species-wide tables",
            "strict_cai": (
                "Requires an a priori, documented set of highly expressed genes "
                "from the relevant strain, cell line or tissue"
            ),
            "expression_yield_predicted": False,
        },
        "hosts": records,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
