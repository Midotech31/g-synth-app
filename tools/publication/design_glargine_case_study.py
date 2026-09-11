#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from gsynth_engine.cloning import clone, translate
from gsynth_engine.constants import ALL_ENZYMES
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.vectors import DEFAULT_VECTOR, sequence_of

CHAINS = {
    "A": {
        "insert": "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA",
        "peptide": "GIVEQCCTSICSLYQLENYCG",
        "recorded_forward": (
            "TATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCT"
            "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAAC"
        ),
        "recorded_reverse": (
            "TCGAGTTAGCCGCAGTAGTTTTCCAGCTGGTACAGGCTGCAGATGCTGGTGCAGCACTGTTCCACGATGCC"
            "AGAACCACGCGGCACCAGACCAGAAGAGTGGTGGTGGTGGTGGTGAGAAGAACCCA"
        ),
    },
    "B": {
        "insert": (
            "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
            "TTTTTTTACACCCCAAAAACCCGCCGCTAA"
        ),
        "peptide": "FVNQHLCGSHLVEALYLVCGERGFFYTPKTRR",
        "recorded_forward": (
            "TATGGGTTCTTCTCACCACCACCACCACCACTCTTCTGGTCTGGTGCCGCGTGGTTCT"
            "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
            "TTTTTTTACACCCCAAAAACCCGCCGCTAAC"
        ),
        "recorded_reverse": (
            "TCGAGTTAGCGGCGGGTTTTTGGGGTGTAAAAAAAGCCGCGTTCGCCGCACACCAGGTACAGCGCTTCC"
            "ACCAGATGGCTGCCGCACAGATGCTGGTTCACAAAAGAACCACGCGGCACCAGACCAGAAGAGTGGTGG"
            "TGGTGGTGGTGAGAAGAACCCA"
        ),
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-record", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    vector_record = sequence_of(DEFAULT_VECTOR.key)
    if vector_record is None:
        raise RuntimeError("The verified pET-21a(+) reference is unavailable.")

    evidence = {
        "analysis": {
            "software": "G-Synth",
            "version": "1.0.0",
            "workflow": "Small Sequence Design followed by in-silico NdeI/XhoI cloning",
            "enzyme_pair": "NdeI / XhoI",
            "cassette": "MGSSHHHHHHSSGLVPRGS",
            "vector": DEFAULT_VECTOR.name,
        },
        "experimental_source_record": {
            "path": str(args.source_record),
            "sha256": sha256(args.source_record),
        },
        "chains": {},
    }

    for name, source in CHAINS.items():
        design = design_small_sequence(
            source["insert"],
            enzyme_pair="NdeI / XhoI",
            is_coding=False,
            cleavage_site="Thrombin",
        )
        if design.forward != source["recorded_forward"]:
            raise AssertionError(f"{name}: G-Synth forward oligo differs from the experimental record")
        if design.reverse != source["recorded_reverse"]:
            raise AssertionError(f"{name}: G-Synth reverse oligo differs from the experimental record")

        bottom_in_top_sense = reverse_complement(design.reverse)
        ndei = ALL_ENZYMES["NdeI"]
        offset = int(ndei["cut_bottom"]) - int(ndei["cut_top"])
        core_top = design.forward[offset:]
        core_bottom = bottom_in_top_sense[:len(core_top)]
        if core_top != core_bottom:
            raise AssertionError(f"{name}: the designed strands do not anneal perfectly")

        coding_dna = design.coding_region
        protein = translate(coding_dna)
        expected_protein = "MGSSHHHHHHSSGLVPRGS" + source["peptide"] + "*"
        if protein != expected_protein:
            raise AssertionError(f"{name}: unexpected translated product {protein}")

        cloned = clone(
            vector_record["sequence"],
            design.forward,
            left_enzyme="NdeI",
            right_enzyme="XhoI",
            name=f"pET-21a(+)-glargine-{name}",
            vector_annotations=vector_record.get("features", []),
            vector_spec=DEFAULT_VECTOR,
            insert_reverse=design.reverse,
            orf_start=design.orf_start,
        )
        if not cloned.is_clonable:
            raise AssertionError(f"{name}: cloning failed: {cloned.problems}")

        evidence["chains"][name] = {
            "input": {
                "optimised_coding_sequence": source["insert"],
                "length_bp": len(source["insert"]),
                "expected_chain_peptide": source["peptide"],
            },
            "design": {
                "forward_oligo_5_to_3": design.forward,
                "reverse_oligo_5_to_3": design.reverse,
                "forward_length_nt": design.forward_length,
                "reverse_length_nt": design.reverse_length,
                "forward_gc_percent": design.forward_gc,
                "reverse_gc_percent": design.reverse_gc,
                "forward_tm_c": design.forward_tm,
                "reverse_tm_c": design.reverse_tm,
                "left_overhang_5_to_3": design.left_overhang,
                "right_overhang_5_to_3": design.right_overhang,
                "orf_start_zero_based": design.orf_start,
                "coding_sequence": coding_dna,
                "translated_product": protein,
                "expected_product_match": protein == expected_protein,
                "experimental_record_forward_exact_match": True,
                "experimental_record_reverse_exact_match": True,
                "duplex_core_exact_match": True,
                "segments": [asdict(segment) for segment in design.segments],
                "warnings": design.warnings,
            },
            "cloning": {
                "is_clonable": cloned.is_clonable,
                "recombinant_length_bp": cloned.length,
                "insert_length_top_strand_nt": cloned.insert_length,
                "insert_start_zero_based": cloned.insert_start,
                "insert_end_zero_based_half_open": cloned.insert_end,
                "reversed_insert_relative_to_vector_numbering": cloned.reversed_insert,
                "junctions": [asdict(junction) for junction in cloned.junctions],
                "translated_product": cloned.protein,
                "problems": cloned.problems,
                "warnings": cloned.warnings,
                "vector_tag_outcomes": [asdict(tag) for tag in cloned.tags],
            },
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
