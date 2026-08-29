"""One preflight contract for PCR, assembly, cloning and verification."""
from __future__ import annotations

from gsynth_engine.diagnostics import PreflightCheck, PreflightReport, make_report
from gsynth_engine.sequence import reverse_complement


def ssd_preflight(result) -> PreflightReport:
    coding_ok = not result.is_coding or result.coding_region.startswith("ATG")
    too_long = max(result.forward_length, result.reverse_length) > 200
    checks = [
        PreflightCheck(
            "SSD_START_CODON",
            "Coding sequence begins at an initiating ATG",
            "pass" if coding_ok else "block",
            f"The coding frame begins at oligo base {result.orf_start + 1}."
            if coding_ok else "The coding frame does not begin with ATG.",
            "Correct the insert start or choose the non-coding cassette mode.",
        ),
        PreflightCheck(
            "SSD_TERMINAL_ENDS",
            "The duplex presents the selected restriction ends",
            "pass",
            f"The insert presents {result.left_enzyme} {result.left_overhang or 'blunt'} "
            f"and {result.right_enzyme} {result.right_overhang or 'blunt'} ends.",
        ),
        PreflightCheck(
            "SSD_OLIGO_LENGTH",
            "Oligos are within the direct-synthesis range",
            "review" if too_long else "pass",
            f"The longer oligo is {max(result.forward_length, result.reverse_length)} nt.",
            "Use fragmented assembly when the supplier cannot synthesize this length.",
        ),
        PreflightCheck(
            "SSD_REVIEW_NOTES",
            "Design notes reviewed",
            "review" if result.warnings else "pass",
            " ".join(result.warnings) or "No additional SSD warnings.",
            "Review internal-site and start-codon notes before ordering.",
        ),
    ]
    return make_report("ssd", checks)


def pcr_preflight(result, *, keep_frame: bool = False) -> PreflightReport:
    cloning = result.left_enzyme is not None and result.right_enzyme is not None
    checks = [
        PreflightCheck(
            "PCR_TARGET_AMPLIFIED",
            "Requested target is represented",
            "pass" if result.amplified_region else "block",
            f"{len(result.amplified_region)} template bases are in the product.",
            "Choose a non-empty target within the template.",
            {"template_start": result.template_start, "template_end": result.template_end},
        ),
    ]
    if cloning:
        checks.extend([
            PreflightCheck(
                "PCR_RESTRICTION_SITE_COUNT",
                "Selected enzymes cut only the intended tails",
                "pass" if not result.problems else "block",
                "Each selected enzyme appears once at its intended product end."
                if not result.problems else " ".join(result.problems),
                "Choose another enzyme pair or remove the internal/overlapping site.",
            ),
            PreflightCheck(
                "PCR_DIGEST_SIMULATED",
                "The PCR product digests into one insert",
                "pass" if result.digest is not None else "block",
                f"A {result.digest.length} nt top strand remains after digestion."
                if result.digest is not None else "No digestible insert was produced.",
                "Resolve the restriction-site check before ordering primers.",
            ),
        ])
    if keep_frame:
        frame_notes = [
            note for note in result.warnings
            if "frame" in note.lower() or "codon" in note.lower()
        ]
        checks.append(PreflightCheck(
            "PCR_READING_FRAME",
            "The requested reading frame is preserved",
            "review" if frame_notes else "pass",
            " ".join(frame_notes) or (
                f"Translation begins at insert base {result.insert_orf_start + 1}."
                if result.insert_orf_start is not None else "Frame preservation was not requested."
            ),
            "Review translated junction bases and any internal stop before ordering.",
        ))
    checks.append(PreflightCheck(
        "PCR_REVIEW_NOTES",
        "Primer and strategy notes reviewed",
        "review" if result.warnings else "pass",
        " ".join(result.warnings) or "No additional primer or strategy warnings.",
        "Confirm informational choices and primer-quality warnings at the bench.",
    ))
    return make_report("pcr", checks)


def optimisation_preflight(result) -> PreflightReport:
    checks = [
        PreflightCheck(
            "OPTIMISE_PROTEIN_INVARIANT",
            "The encoded protein is preserved",
            "block" if any("protein" in problem.lower() for problem in result.problems) else "pass",
            "The optimized gene translates to the same protein."
            if not any("protein" in problem.lower() for problem in result.problems)
            else " ".join(result.problems),
            "Do not use an optimized sequence whose translation differs from the input.",
        ),
        PreflightCheck(
            "OPTIMISE_CONSTRAINTS",
            "Required sequence constraints are satisfied",
            "block" if result.problems else "pass",
            "No blocking motif, GC, repeat or enzyme-site constraint remains."
            if not result.problems else " ".join(result.problems),
            "Relax a conflicting constraint or change the coding sequence strategy.",
        ),
        PreflightCheck(
            "OPTIMISE_REVIEW_NOTES",
            "Optimization trade-offs reviewed",
            "review" if result.warnings else "pass",
            " ".join(result.warnings) or "No additional optimization warnings.",
            "Review residual rare codons and local sequence-quality trade-offs.",
        ),
    ]
    return make_report("codon_optimisation", checks)


def assembly_preflight(plan) -> PreflightReport:
    verification = plan.verify()
    overhangs = list(plan.junction_overhangs)
    reverse_overhangs = {reverse_complement(overhang) for overhang in overhangs}
    orthogonal = len(set(overhangs)) == len(overhangs) and not any(
        overhang in reverse_overhangs - {reverse_complement(overhang)}
        for overhang in overhangs
    )
    checks = [
        PreflightCheck(
            "ASSEMBLY_BOTH_STRANDS_RECONSTRUCT",
            "Both strands reconstruct the designed cassette",
            "pass" if not verification else "block",
            "Re-ligation reproduces the forward and reverse designs base for base."
            if not verification else " ".join(verification),
            "Do not order the oligos until fragment reconstruction is exact.",
        ),
        PreflightCheck(
            "ASSEMBLY_JUNCTIONS_ORTHOGONAL",
            "Internal junction overhangs are unique",
            "pass" if orthogonal else "block",
            f"{len(overhangs)} internal overhangs were checked for reuse and reverse complements.",
            "Redesign the fragment boundaries with orthogonal junctions.",
            {"overhangs": overhangs},
        ),
        PreflightCheck(
            "ASSEMBLY_TERMINAL_ENDS",
            "Terminal ends match the selected restriction strategy",
            "pass",
            f"The assembled insert presents {plan.terminal_ends[0]} and {plan.terminal_ends[1]}.",
        ),
        PreflightCheck(
            "ASSEMBLY_REVIEW_NOTES",
            "Assembly notes reviewed",
            "review" if plan.warnings else "pass",
            " ".join(plan.warnings) or "No additional assembly warnings.",
            "Review synthesis and annealing notes before ordering.",
        ),
    ]
    return make_report("assembly", checks)


def cloning_preflight(result, *, duplex_mismatches: list[int] | None = None) -> PreflightReport:
    mismatches = duplex_mismatches or []
    end_problems = [problem for problem in result.problems if "does not match" in problem]
    duplex_problems = [problem for problem in result.problems if "do not pair" in problem]
    frame_problems = [problem for problem in result.problems if "truncated" in problem]
    regenerated = all(junction.site_regenerated for junction in result.junctions)
    start_ok = not result.protein or result.protein.startswith("M")
    checks = [
        PreflightCheck(
            "CLONE_END_COMPATIBILITY",
            "Insert and vector ends are compatible",
            "pass" if not end_problems else "block",
            "Both observed insert ends anneal to the cut vector."
            if not end_problems else " ".join(end_problems),
            "Use an insert digested with the selected enzymes or redesign its ends.",
        ),
        PreflightCheck(
            "CLONE_DUPLEX_INTEGRITY",
            "Both insert strands pair everywhere",
            "pass" if not mismatches and not duplex_problems else "block",
            "No mismatch in the supplied insert duplex."
            if not mismatches and not duplex_problems
            else " ".join(duplex_problems) or f"{len(mismatches)} insert positions do not pair.",
            "Correct the oligo strands before ligation.",
            {"mismatch_positions": mismatches},
        ),
        PreflightCheck(
            "CLONE_VECTOR_SITE_COUNT",
            "Each enzyme cuts the vector once",
            "pass",
            "The vector linearisation gate confirmed one site per enzyme.",
        ),
        PreflightCheck(
            "CLONE_ORIENTATION_FORCED",
            "Insert orientation is forced",
            "pass" if result.left_enzyme != result.right_enzyme else "block",
            f"Different {result.left_enzyme} and {result.right_enzyme} ends define orientation.",
            "Use two different restriction ends.",
        ),
        PreflightCheck(
            "CLONE_JUNCTION_SITES",
            "Restriction sites regenerate at both junctions",
            "pass" if regenerated else "review",
            "Both sites regenerate, so the insert can be excised diagnostically."
            if regenerated else "At least one cloning site is lost after ligation.",
            "Choose whether loss of the diagnostic site is acceptable.",
        ),
        PreflightCheck(
            "CLONE_READING_FRAME",
            "Translation starts correctly and reaches the intended terminus",
            "block" if frame_problems else "review" if not start_ok else "pass",
            " ".join(frame_problems) if frame_problems else (
                f"{len(result.protein)} residues translate from an initiating methionine."
                if result.protein else "No reading frame was supplied for translation."
            ),
            "Inspect the ATG, retained junction bases and stop-codon placement.",
        ),
        PreflightCheck(
            "CLONE_REVIEW_NOTES",
            "Cloning and expression notes reviewed",
            "review" if result.warnings else "pass",
            " ".join(result.warnings) or "No additional cloning warnings.",
            "Review diagnostic-digest, tag and stop-codon notes before proceeding.",
        ),
    ]
    return make_report("cloning", checks)


def verification_state(report) -> str:
    """Stable five-state sequencing verdict used by API and interface."""
    if report.is_verified:
        return "fully_verified"
    if report.differences:
        return "differences_detected"
    if report.reads:
        return "partial_match"
    if report.warnings:
        return "reads_unplaced"
    return "not_checked"


def verification_preflight(report) -> PreflightReport:
    checks = [
        PreflightCheck(
            "VERIFY_READS_PLACED",
            "At least one sequencing read was placed",
            "pass" if report.reads else "block",
            f"{len(report.reads)} reads aligned to the design."
            if report.reads else "No sequencing read could be aligned.",
            "Check read orientation, trim settings and that the intended construct was sequenced.",
        ),
        PreflightCheck(
            "VERIFY_COMPLETE_COVERAGE",
            "The requested region is fully covered",
            "pass" if report.fully_covered else "block",
            f"{report.coverage}% of the requested region is covered.",
            "Add sequencing reads across every remaining gap.",
            {"gaps": list(report.gaps), "coverage": report.coverage},
        ),
        PreflightCheck(
            "VERIFY_SEQUENCE_AGREEMENT",
            "All covered bases agree with the design",
            "pass" if not report.differences else "block",
            "No differences were observed."
            if not report.differences else f"{len(report.differences)} differences were observed.",
            "Inspect trace quality and confirm or reject each sequence difference.",
        ),
    ]
    result = make_report("verification", checks)
    # State is a workflow fact rather than a sixth verdict vocabulary.
    return result
