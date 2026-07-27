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
from gsynth_engine.merzoug import AssemblyPlan, design_merzoug_assembly
from gsynth_engine.protocol import bench_protocol, order_sheet, order_sheet_csv
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import SSDResult, design_small_sequence
from rest_framework import status
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.design.serializers import (
    CLEAVAGE_NAMES,
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
                "is_first": fragment.is_first,
                "is_last": fragment.is_last,
            }
            for fragment in plan.fragments
        ],
        "oligos": [order.as_row for order in order_sheet(plan, construct_name=construct_name)],
        "ssd": _ssd_payload(plan.ssd),
        "warnings": plan.warnings,
        # Empty means: re-ligating these oligos reproduces the design exactly.
        "verification": plan.verify(),
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
