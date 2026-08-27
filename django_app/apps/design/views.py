"""Design endpoints — a thin HTTP layer over `gsynth_engine`.

Every view here follows the same shape: validate the request, call the
engine, translate a `SequenceError` into a 400 the user can act on, and
serialise the result. The biology stays in the engine, which is where the
tests are.
"""
from __future__ import annotations

from django.http import HttpResponse
from rest_framework import status
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.design.serializers import (
    CLEAVAGE_NAMES,
    AlignRequestSerializer,
    CloneRequestSerializer,
    LigationRequestSerializer,
    OptimiseRequestSerializer,
    PcrRequestSerializer,
    PrimerRequestSerializer,
    SaveableAssemblyRequestSerializer,
    SSDRequestSerializer,
    TraceUploadSerializer,
    VerifyRequestSerializer,
    resolve_vector,
)
from apps.projects.models import Project
from gsynth_engine import vectors as vector_catalogue
from gsynth_engine.align import Scoring, align, blosum62
from gsynth_engine.chromatogram import read_ab1, summarise
from gsynth_engine.cloning import (
    CloningResult,
    clone,
    find_sites,
    open_reading_frames,
)
from gsynth_engine.codon import (
    ECOLI,
    Constraints,
    OptimisationResult,
    build_table,
    optimise,
)
from gsynth_engine.constants import (
    ALL_ENZYMES,
    CLEAVAGE_SITES,
    COMMON_ENZYME_PAIRS,
    RESTRICTION_ENZYMES,
    overhang,
    supplies_start_codon,
)
from gsynth_engine.duplex import DuplexView, construct_duplex, junction_view
from gsynth_engine.genbank import oligos_to_fasta, to_fasta, to_genbank
from gsynth_engine.ligation import ligation_series, plan_ligation
from gsynth_engine.merzoug import AssemblyPlan, design_merzoug_assembly
from gsynth_engine.pcr import design_pcr
from gsynth_engine.preflight import (
    assembly_preflight,
    cloning_preflight,
    optimisation_preflight,
    pcr_preflight,
    ssd_preflight,
    verification_preflight,
    verification_state,
)
from gsynth_engine.primers import design_sequencing_primers
from gsynth_engine.protocol import (
    bench_protocol,
    cloning_worksheet,
    order_sheet,
    order_sheet_csv,
)
from gsynth_engine.provenance import build_provenance
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import SSDResult, design_small_sequence
from gsynth_engine.thermo import ANNEALING
from gsynth_engine.verify import verify


def _bad_request(error: SequenceError) -> Response:
    """Engine errors are already written for the user — pass them through."""
    return Response({"detail": str(error)}, status=status.HTTP_400_BAD_REQUEST)


def _provenance(workflow: str, data: dict, output: str, *, vector: str = "") -> dict:
    """Timestamped reproducibility record for responses, saves and exports."""
    from django.utils import timezone

    return build_provenance(
        workflow,
        parameters=dict(data),
        output_sequence=output,
        vector_sequence=vector,
        generated_at=timezone.now().isoformat(),
    )


def _ssd_payload(result: SSDResult) -> dict:
    return {
        "forward": result.forward,
        "reverse": result.reverse,
        "forward_length": result.forward_length,
        "reverse_length": result.reverse_length,
        "forward_gc": result.forward_gc,
        "reverse_gc": result.reverse_gc,
        "forward_tm": result.forward_tm,
        "reverse_tm": result.reverse_tm,
        "left_enzyme": result.left_enzyme,
        "right_enzyme": result.right_enzyme,
        "left_overhang": result.left_overhang,
        "right_overhang": result.right_overhang,
        "cleavage_site": result.cleavage_site,
        "orf_start": result.orf_start,
        "coding_region": result.coding_region,
        "segments": [
            {
                "name": segment.name,
                "start": segment.start,
                "end": segment.end,
                "sequence": segment.sequence,
            }
            for segment in result.segments
        ],
        "warnings": result.warnings,
        "preflight": ssd_preflight(result).to_dict(),
    }


def _duplex_payload(view: DuplexView) -> dict:
    """The hybridisation view, as coordinates the client draws from.

    Sent as two padded strings plus spans rather than as pre-wrapped lines,
    so the browser can re-wrap to its own width without asking again.
    """
    return {
        "top": view.top,
        "bottom": view.bottom,
        "pairs": view.paired(),
        "width": view.width,
        "left_overhang": view.left_overhang,
        "right_overhang": view.right_overhang,
        "junctions": view.junctions,
        "mismatches": view.mismatches(),
        "segments": [
            {"name": span.name, "start": span.start, "end": span.end}
            for span in view.segments
        ],
        "top_fragments": [
            {"name": span.name, "start": span.start, "end": span.end}
            for span in view.top_fragments
        ],
        "bottom_fragments": [
            {"name": span.name, "start": span.start, "end": span.end}
            for span in view.bottom_fragments
        ],
    }


def _assembly_payload(plan: AssemblyPlan, construct_name: str) -> dict:
    return {
        "construct_forward": plan.construct_forward,
        "construct_reverse": plan.construct_reverse,
        "construct_length": plan.construct_length,
        "construct_gc": round(gc_content(plan.construct_forward), 1),
        "fragment_count": plan.fragment_count,
        "oligo_count": plan.oligo_count,
        "overhang_length": plan.overhang_length,
        "longest_oligo": plan.longest_oligo,
        "junction_overhangs": plan.junction_overhangs,
        # Measured off the assembled fragments, not copied from the design:
        # this is what the vector will actually be offered.
        "terminal_ends": [
            {"side": side, "enzyme": enzyme, "overhang": sequence, "kind": kind}
            for side, enzyme, (sequence, kind) in (
                ("left", plan.ssd.left_enzyme, plan.terminal_ends[0]),
                ("right", plan.ssd.right_enzyme, plan.terminal_ends[1]),
            )
        ],
        "fragments": [
            {
                "index": fragment.index,
                "name": fragment.name,
                "forward": fragment.forward,
                "reverse": fragment.reverse,
                "forward_length": fragment.forward_length,
                "reverse_length": fragment.reverse_length,
                "forward_tm": fragment.forward_tm,
                "reverse_tm": fragment.reverse_tm,
                "top_start": fragment.top_start,
                "top_end": fragment.top_end,
                "left_overhang": fragment.left_overhang,
                "right_overhang": fragment.right_overhang,
                "left_overhang_strand": fragment.left_overhang_strand,
                "right_overhang_strand": fragment.right_overhang_strand,
                "bottom_offset": fragment.bottom_offset,
                "is_first": fragment.is_first,
                "is_last": fragment.is_last,
            }
            for fragment in plan.fragments
        ],
        "oligos": [order.as_row for order in order_sheet(plan, construct_name=construct_name)],
        "ssd": _ssd_payload(plan.ssd),
        "duplex": _duplex_payload(construct_duplex(plan)),
        # Tm is meaningless without the reaction it refers to.
        "tm_conditions": {
            "name": ANNEALING.name,
            "summary": ANNEALING.summary,
            "model": "Nearest-neighbour (SantaLucia 1998)",
        },
        "warnings": plan.warnings,
        # Empty means: re-ligating these oligos reproduces the design exactly.
        "verification": plan.verify(),
        "preflight": assembly_preflight(plan).to_dict(),
    }


#: Enzymes worth marking on a recombinant map — every one with verified cut
#: geometry, because "does anything else cut here" is the question that decides
#: whether a diagnostic digest will work.
def _restriction_annotations(plasmid: str, highlight: tuple[str, ...]) -> list[dict]:
    """Every single-cutter, plus the pair used, as drawable features.

    Single cutters only: a site that appears eleven times is noise on a map
    and useless for a diagnostic digest. The two cloning enzymes are marked
    even when they cut more than once, because that is exactly the case the
    user needs to see.
    """
    out: list[dict] = []
    for enzyme in sorted(ALL_ENZYMES):
        sites = find_sites(plasmid, enzyme, circular=True)
        if not sites:
            continue
        used = enzyme in highlight
        if len(sites) > 1 and not used:
            continue
        site = str(ALL_ENZYMES[enzyme]["recognition"])
        for position in sites:
            end = position + len(site)
            out.append({
                "name": enzyme,
                "type": "restriction_site",
                "start": position,
                "end": end,
                "direction": 1,
                "color": "#9E3D3D" if used else "#78889B",
                "cuts": len(sites),
                "used": used,
                "recognition": site,
                # A circular molecule has no beginning, so a site can straddle
                # position 0. The coordinates stay honest and say so, rather
                # than being clamped into a feature that is not the site.
                "wraps": end > len(plasmid),
            })
    return out


def _validation(result: CloningResult, duplex_mismatches: list[int]) -> list[dict]:
    """The checks, each stated as a claim that passed or did not.

    One banner collapses a dozen distinct questions into a colour. Listing
    them separately means a design that fails on orientation and a design
    that fails on a duplicated site do not look alike.
    """
    junctions_ok = all(j.site_regenerated for j in result.junctions)
    ends_ok = not any("does not match" in p for p in result.problems)
    duplex_problems = [problem for problem in result.problems if "do not pair" in problem]
    frame_ok = not any("truncated" in p for p in result.problems)

    return [
        {
            "check": "Overhangs are compatible",
            "passed": ends_ok,
            "detail": " ".join(p for p in result.problems if "does not match" in p)
            or f"{result.left_enzyme} and {result.right_enzyme} ends match the vector.",
        },
        {
            "check": "Both strands pair everywhere",
            "passed": not duplex_mismatches and not duplex_problems,
            "detail": (
                f"{len(duplex_mismatches)} positions do not pair."
                if duplex_mismatches else " ".join(duplex_problems)
                if duplex_problems else "No mismatch in the insert duplex."
            ),
        },
        {
            "check": "Each enzyme cuts the vector once",
            "passed": True,          # a second site raises before we get here
            "detail": "Checked when the vector was cut; a second site is refused.",
        },
        {
            "check": "Orientation is forced",
            "passed": result.left_enzyme != result.right_enzyme,
            "detail": (
                "The insert reads on the minus strand of the vector's own "
                "numbering, as every pET cassette does."
                if result.reversed_insert
                else "The insert reads with the vector's numbering."
            ),
        },
        {
            "check": "Sites are regenerated at both seams",
            "passed": junctions_ok,
            "detail": (
                "The insert can be cut back out — which is how the clone gets "
                "verified on a gel."
                if junctions_ok else
                "At least one site is lost, so the insert cannot be excised."
            ),
        },
        {
            "check": "Reading frame survives the junction",
            "passed": frame_ok,
            "detail": " ".join(p for p in result.problems if "truncated" in p)
            or (f"{len(result.protein)} residues translated."
                if result.protein else "No reading frame was given."),
        },
    ]


def _clone_payload(result: CloningResult, ssd: SSDResult | None, plan) -> dict:
    """The recombinant plasmid, shaped for a map viewer.

    The insert is sent as one more annotation so the client can draw the
    whole plasmid from a single list, rather than special-casing it.
    """
    annotations = list(result.annotations)
    annotations.append({
        "name": result.name,
        "type": "CDS" if result.protein else "misc_feature",
        "start": result.insert_start,
        "end": result.insert_end,
        "direction": 1,
        "color": "#0E6E77",
    })

    duplex_mismatches = construct_duplex(plan).mismatches() if plan else []
    return {
        "plasmid": result.plasmid,
        "name": result.name,
        "length": result.length,
        "gc": result.gc,
        "topology": "circular",
        "insert_start": result.insert_start,
        "insert_end": result.insert_end,
        "insert_length": result.insert_length,
        "backbone_length": result.backbone_length,
        "removed_length": result.removed_length,
        "left_enzyme": result.left_enzyme,
        "right_enzyme": result.right_enzyme,
        "protein": result.protein,
        "protein_length": len(result.protein),
        # Every pET cassette reads on the minus strand of the supplier's
        # numbering, so this is the normal case rather than a warning.
        "reversed_insert": result.reversed_insert,
        "tags": [
            {
                "name": tag.name,
                "end": tag.end,
                "present": tag.present,
                "position": tag.position,
                "note": tag.note,
            }
            for tag in result.tags
        ],
        "annotations": annotations,
        "junctions": [
            {
                "name": junction.name,
                "enzyme": junction.enzyme,
                "overhang": junction.overhang,
                "kind": junction.kind,
                "position": junction.position,
                "context": junction.context,
                "site_regenerated": junction.site_regenerated,
            }
            for junction in result.junctions
        ],
        "orfs": open_reading_frames(result.plasmid, minimum_codons=40)[:5],
        # Each seam drawn as the two ends that made it, so "the overhangs
        # match" can be checked rather than believed.
        "junction_views": [
            {
                "name": view.name,
                "enzyme": view.enzyme,
                "overhang": view.overhang,
                "kind": view.kind,
                "compatible": view.compatible,
                "reason": view.reason,
                "left_top": view.left_top,
                "left_bottom": view.left_bottom,
                "right_top": view.right_top,
                "right_bottom": view.right_bottom,
                "joined_top": view.joined_top,
                "joined_bottom": view.joined_bottom,
                "joined_pairs": view.joined_pairs,
                "seam": view.seam,
                "overhang_span": list(view.overhang_span),
            }
            for view in (
                junction_view(
                    result.plasmid, name=j.name, enzyme=j.enzyme,
                    position=j.position, overhang=j.overhang, kind=j.kind,
                    strand="", flank=18,
                )
                for j in result.junctions
            )
        ],
        "restriction_sites": _restriction_annotations(
            result.plasmid, (result.left_enzyme, result.right_enzyme),
        ),
        "validation": _validation(
            result,
            duplex_mismatches,
        ),
        "preflight": cloning_preflight(
            result, duplex_mismatches=duplex_mismatches,
        ).to_dict(),
        "warnings": result.warnings,
        # Empty means these two molecules really do join.
        "problems": result.problems,
        "is_clonable": result.is_clonable,
        # None when the caller supplied an insert that was already cut: there
        # was no SSD design, and inventing