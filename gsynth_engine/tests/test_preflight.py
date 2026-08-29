from gsynth_engine import vectors
from gsynth_engine.cloning import clone
from gsynth_engine.codon import optimise
from gsynth_engine.merzoug import design_merzoug_assembly
from gsynth_engine.pcr import design_pcr
from gsynth_engine.preflight import (
    assembly_preflight,
    cloning_preflight,
    optimisation_preflight,
    pcr_preflight,
    verification_preflight,
    verification_state,
)
from gsynth_engine.verify import verify

GENE = (
    "ATGAAAGGTGAAGAATTGTTCACCGGTGTTGTTCCGATTCTGGTTGAACTGGATGGTGATGTT"
    "AACGGTCACAAATTCTCTGTTTCTGGTGAAGGTGAAGGTGATGCTACCTACGGTAAACTGACC"
    "CTGAAATAA"
)


def test_every_workflow_uses_the_same_preflight_shape():
    pcr = design_pcr(GENE, left_enzyme="NdeI", right_enzyme="XhoI", keep_frame=True)
    assembly = design_merzoug_assembly(GENE, is_coding=True)
    plasmid = clone(
        vectors.sequence_of("pET-21a")["sequence"], assembly.construct_forward,
        insert_reverse=assembly.construct_reverse,
        left_enzyme="NdeI", right_enzyme="XhoI", orf_start=assembly.ssd.orf_start,
    )
    reports = [
        pcr_preflight(pcr, keep_frame=True),
        assembly_preflight(assembly),
        cloning_preflight(plasmid),
    ]
    for report in reports:
        payload = report.to_dict()
        assert payload["verdict"] in {"ready", "review", "blocked"}
        assert isinstance(payload["can_export"], bool)
        assert all({"code", "label", "status", "detail", "passed"} <= check.keys()
                   for check in payload["checks"])


def test_a_broken_pcr_is_blocked_with_a_stable_code():
    gene = "ATG" + "GCTAGC" + GENE[3:]
    report = pcr_preflight(
        design_pcr(gene, left_enzyme="NheI", right_enzyme="XhoI")
    )
    assert report.verdict == "blocked"
    assert "PCR_RESTRICTION_SITE_COUNT" in {
        diagnostic.code for diagnostic in report.diagnostics
    }


def test_optimisation_uses_the_same_release_contract():
    result = optimise(GENE)
    report = optimisation_preflight(result)
    assert report.workflow == "codon_optimisation"
    assert "OPTIMISE_PROTEIN_INVARIANT" in {check.code for check in report.checks}
    assert report.verdict != "blocked"


def test_verification_has_five_explicit_states():
    partial = verify(GENE, {"read": GENE[:70]}, trim=0)
    complete = verify(GENE, {"read": GENE}, trim=0)
    absent = verify(GENE, {}, trim=0)
    unplaced = verify(GENE, {"junk": "ACGT" * 30}, trim=0)
    changed = list(GENE)
    changed[40] = "A" if changed[40] != "A" else "C"
    different = verify(GENE, {"read": "".join(changed)}, trim=0)

    assert verification_state(complete) == "fully_verified"
    assert verification_state(partial) == "partial_match"
    assert verification_state(absent) == "not_checked"
    assert verification_state(unplaced) == "reads_unplaced"
    assert verification_state(different) == "differences_detected"
    assert verification_preflight(partial).verdict == "blocked"
