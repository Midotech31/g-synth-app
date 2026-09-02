"""Generate machine-readable evidence for peptide starts and enzyme coverage."""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from gsynth_engine.cloning import clone, translate
from gsynth_engine.codon import TABLES, Constraints, optimise
from gsynth_engine.constants import ALL_ENZYMES, overhang
from gsynth_engine.sequence import reverse_complement
from gsynth_engine.ssd import design_small_sequence
from gsynth_engine.vectors import DEFAULT_VECTOR, sequence_of

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "publication_evidence" / "peptide_and_enzyme_validation.json"
MATURE_GLARGINE_A = "GIVEQCCTSICSLYQLENYCG"
MATURE_GLARGINE_B = "FVNQHLCGSHLVEALYLVCGERGFFYTPKTRR"
EXPERIMENTAL_CODING_SEQUENCES = {
    "A": "GGCATCGTGGAACAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAAAACTACTGCGGCTAA",
    "B": (
        "TTTGTGAACCAGCATCTGTGCGGCAGCCATCTGGTGGAAGCGCTGTACCTGGTGTGCGGCGAACGCGGC"
        "TTTTTTTACACCCCAAAAACCCGCCGCTAA"
    ),
}


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


def de_novo_glargine_design(chain: str, peptide: str) -> dict:
    constraints = Constraints(
        avoid_enzymes=("NdeI", "XhoI"),
        max_homopolymer=5,
        gc_min=40,
        gc_max=60,
        gc_window=50,
        max_repeat=15,
        avoid_rare=True,
    )
    optimised = optimise(
        peptide,
        table=TABLES["ecoli"],
        constraints=constraints,
        is_protein=True,
        protein_context="auto",
        keep_stop=True,
    )
    design = design_small_sequence(
        optimised.sequence,
        enzyme_pair="NdeI / XhoI",
        is_coding=optimised.recommended_design_is_coding,
        cleavage_site="Thrombin",
    )

    bottom_in_top_sense = reverse_complement(design.reverse)
    ndei = ALL_ENZYMES["NdeI"]
    left_displacement = int(ndei["cut_bottom"]) - int(ndei["cut_top"])
    paired_top = design.forward[left_displacement:]
    paired_bottom = bottom_in_top_sense[:len(paired_top)]
    duplex_core_exact = paired_top == paired_bottom

    vector_record = sequence_of(DEFAULT_VECTOR.key)
    if vector_record is None:
        raise RuntimeError("The verified pET-21a(+) reference is unavailable.")
    cloned = clone(
        vector_record["sequence"],
        design.forward,
        left_enzyme="NdeI",
        right_enzyme="XhoI",
        name=f"pET-21a(+)-de-novo-glargine-{chain}",
        vector_annotations=vector_record.get("features", []),
        vector_spec=DEFAULT_VECTOR,
        insert_reverse=design.reverse,
        orf_start=design.orf_start,
    )

    experimental_sequence = EXPERIMENTAL_CODING_SEQUENCES[chain]
    translated_cassette = translate(design.coding_region)
    expected_cassette = "MGSSHHHHHHSSGLVPRGS" + peptide + "*"
    if not all((
        optimised.is_clean,
        optimised.protein_context == "mature_peptide",
        not optimised.initiator_methionine_added,
        translate(optimised.sequence) == peptide + "*",
        duplex_core_exact,
        translated_cassette == expected_cassette,
        cloned.is_clonable,
    )):
        raise AssertionError(f"{chain}: de-novo peptide-to-clone invariant failed")

    return {
        "input_peptide": peptide,
        "requested_peptide_role": "auto",
        "resolved_peptide_role": optimised.protein_context,
        "initiator_methionine_added": optimised.initiator_methionine_added,
        "expression_host": {
            "key": "ecoli",
            "name": TABLES["ecoli"].name,
            "source": TABLES["ecoli"].source,
        },
        "back_translated_coding_sequence_5_to_3": optimised.sequence,
        "coding_length_bp": optimised.length,
        "profile_relative_cai": optimised.cai_after,
        "gc_percent": optimised.gc_after,
        "low_frequency_codons": optimised.rare_codons_after,
        "translation_preserved": translate(optimised.sequence) == peptide + "*",
        "internal_nde_i_or_xho_i_absent": all(
            ALL_ENZYMES[name]["recognition"] not in optimised.sequence
            for name in ("NdeI", "XhoI")
        ),
        "sequence_sha256": digest(optimised.sequence),
        "experimental_sequence_sha256": digest(experimental_sequence),
        "synonymous_sequence_differs_from_experimental": optimised.sequence != experimental_sequence,
        "experimental_translation_matches_peptide": translate(experimental_sequence) == peptide + "*",
        "ssd": {
            "forward_order_molecule_5_to_3": design.forward,
            "reverse_order_molecule_5_to_3": design.reverse,
            "forward_length_nt": design.forward_length,
            "reverse_length_nt": design.reverse_length,
            "duplex_core_exact": duplex_core_exact,
            "left_overhang_5_to_3": design.left_overhang,
            "right_overhang_5_to_3": design.right_overhang,
            "translated_cassette": translated_cassette,
            "segments": [asdict(segment) for segment in design.segments],
        },
        "restriction_cloning": {
            "vector": DEFAULT_VECTOR.name,
            "is_clonable": cloned.is_clonable,
            "recombinant_length_bp": cloned.length,
            "junctions": [asdict(junction) for junction in cloned.junctions],
            "translated_product": cloned.protein,
        },
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
        "schema_version": "1.1",
        "de_novo_glargine_workflow": {
            "interpretation": (
                "Current G-Synth demonstration from amino-acid input; the synonymous DNA "
                "sequences differ from the physically synthesized experimental constructs and "
                "are therefore validated computationally rather than experimentally."
            ),
            "A": de_novo_glargine_design("A", MATURE_GLARGINE_A),
            "B": de_novo_glargine_design("B", MATURE_GLARGINE_B),
        },
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
