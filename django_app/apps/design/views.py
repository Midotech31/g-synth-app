"""Design endpoints — a thin HTTP layer over `gsynth_engine`.

Every view here follows the same shape: validate the request, call the
engine, translate a `SequenceError` into a 400 the user can act on, and
serialise the result. The biology stays in the engine, which is where the
tests are.
"""
from __future__ import annotations

from django.http import HttpResponse
from gsynth_engine.constants import (
    CLEAVAGE_SITES,
    COMMON_ENZYME_PAIRS,
    RESTRICTION_ENZYMES,
    overhang,
)
from gsynth_engine.cloning import CloningResult, clone, open_reading_frames
from gsynth_engine.duplex import DuplexView, construct_duplex
from gsynth_engine.merzoug import AssemblyPlan, design_merzoug_assembly
from gsynth_engine.protocol import bench_protocol, order_sheet, order_sheet_csv
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import SSDResult, design_small_sequence
from gsynth_engine.thermo import ANNEALING
from rest_framework import status
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.design.serializers import (
    CLEAVAGE_NAMES,
    CloneRequestSerializer,
    SaveableAssemblyRequestSerializer,
    SSDRequestSerializer,
)
from apps.projects.models import Project


def _bad_request(error: SequenceError) -> Response:
    """Engine errors are already written for the user — pass them through."""
    return Response({"detail": str(error)}, status=status.HTTP_400_BAD_REQUEST)


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
    }


def _clone_payload(result: CloningResult, ssd: SSDResult, plan) -> dict:
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
        "warnings": result.warnings,
        # Empty means these two molecules really do join.
        "problems": result.problems,
        "is_clonable": result.is_clonable,
        "insert": _ssd_payload(ssd),
        "assembly": _assembly_payload(plan, result.name) if plan else None,
    }


class EnzymeCatalogueView(APIView):
    """GET /api/design/enzymes/ — what the UI needs to build its dropdowns.

    Public: it is a reference table, and the sign-up screen may want to show
    it before anyone has an account.
    """

    permission_classes = (AllowAny,)

    def get(self, request):
        enzymes = []
        for name in sorted(RESTRICTION_ENZYMES):
            sequence, kind = overhang(name)
            enzymes.append({
                "name": name,
                "recognition": RESTRICTION_ENZYMES[name]["recognition"],
                "overhang": sequence,
                "overhang_type": kind,
                # NdeI's site contains the ATG, which changes how the insert
                # is built — the UI should say so next to the option.
                "supplies_start_codon": "ATG" in str(RESTRICTION_ENZYMES[name]["recognition"]),
            })
        return Response({
            "enzymes": enzymes,
            "common_pairs": list(COMMON_ENZYME_PAIRS),
            "cleavage_sites": [
                {"name": name, "sequence": CLEAVAGE_SITES[name]} for name in CLEAVAGE_NAMES
            ],
        })


class SSDDesignView(APIView):
    """POST /api/design/ssd/ — one insert, two oligos."""

    def post(self, request):
        serializer = SSDRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        try:
            result = design_small_sequence(
                serializer.validated_data["sequence"], **serializer.engine_kwargs
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _ssd_payload(result)

        if serializer.validated_data["save_as_project"]:
            project = Project.objects.create(
                user=request.user,
                name=serializer.validated_data["name"],
                module="ssd",
                sequence=result.forward,
                notes=f"{result.left_enzyme}/{result.right_enzyme}"
                      + (f" · {result.cleavage_site}" if result.cleavage_site else ""),
                data=payload,
            )
            payload["project_id"] = project.id
        return Response(payload)


class MerzougAssemblyView(APIView):
    """POST /api/design/assembly/ — one insert, an ordered set of oligo pairs."""

    def post(self, request):
        serializer = SaveableAssemblyRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        name = serializer.validated_data["name"]
        try:
            plan = design_merzoug_assembly(
                serializer.validated_data["sequence"], **serializer.engine_kwargs
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _assembly_payload(plan, name)

        if serializer.validated_data["save_as_project"]:
            project = Project.objects.create(
                user=request.user,
                name=name,
                module="merzoug_assembly",
                sequence=plan.construct_forward,
                notes=f"{plan.fragment_count} fragments · "
                      f"{plan.oligo_count} oligos · "
                      f"{plan.overhang_length} nt overhangs",
                data=payload,
            )
            payload["project_id"] = project.id
        return Response(payload)


class CloneView(APIView):
    """POST /api/design/clone/ — design an insert and put it in a vector.

    Returns the recombinant plasmid: sequence, junctions, the protein that
    will be expressed, and the vector's annotations at their new coordinates.
    A design that cannot be cloned comes back with `problems` filled in and
    HTTP 200 — the user needs to see what does not fit, not an error page.
    """

    def post(self, request):
        serializer = CloneRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        name = data["name"]

        try:
            if data["fragment"]:
                plan = design_merzoug_assembly(
                    data["sequence"], **serializer.engine_kwargs
                )
                ssd = plan.ssd
                insert_forward = plan.construct_forward
                insert_reverse = plan.construct_reverse
            else:
                plan = None
                ssd = design_small_sequence(
                    data["sequence"],
                    **{k: v for k, v in serializer.engine_kwargs.items()
                       if k not in ("target_oligo_length", "overhang_length")},
                )
                insert_forward = ssd.forward
                insert_reverse = ssd.reverse

            result = clone(
                data["vector"],
                insert_forward,
                insert_reverse=insert_reverse,
                left_enzyme=data["left_enzyme"],
                right_enzyme=data["right_enzyme"],
                circular=data["vector_is_circular"],
                name=name,
                vector_annotations=[
                    dict(feature) for feature in data.get("vector_annotations") or []
                ],
                orf_start=ssd.orf_start,
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _clone_payload(result, ssd, plan)
        payload["vector_name"] = data["vector_name"]

        if data["save_as_project"] and result.is_clonable:
            project = Project.objects.create(
                user=request.user,
                name=name,
                module="cloning",
                sequence=result.plasmid,
                notes=f"{result.length} bp in {data['vector_name']} · "
                      f"{result.left_enzyme}/{result.right_enzyme}",
                data=payload,
            )
            payload["project_id"] = project.id
        return Response(payload)


class OrderSheetView(APIView):
    """POST /api/design/assembly/order-sheet/ — the oligo list as CSV."""

    def post(self, request):
        serializer = SaveableAssemblyRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        name = serializer.validated_data["name"]
        try:
            plan = design_merzoug_assembly(
                serializer.validated_data["sequence"], **serializer.engine_kwargs
            )
        except SequenceError as error:
            return _bad_request(error)

        csv_text = order_sheet_csv(plan, construct_name=name)
        response = HttpResponse(csv_text, content_type="text/csv")
        safe_name = name.replace(" ", "_") or "construct"
        response["Content-Disposition"] = f'attachment; filename="{safe_name}_oligos.csv"'
        return response


class ProtocolView(APIView):
    """POST /api/design/assembly/protocol/ — the bench protocol as text."""

    def post(self, request):
        serializer = SaveableAssemblyRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        name = serializer.validated_data["name"]
        try:
            plan = design_merzoug_assembly(
                serializer.validated_data["sequence"], **serializer.engine_kwargs
            )
        except SequenceError as error:
            return _bad_request(error)

        text = bench_protocol(plan, construct_name=name)
        response = HttpResponse(text, content_type="text/plain; charset=utf-8")
        safe_name = name.replace(" ", "_") or "construct"
        response["Content-Disposition"] = f'attachment; filename="{safe_name}_protocol.txt"'
        return response
