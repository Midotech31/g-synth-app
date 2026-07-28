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
from gsynth_engine.codon import (
    ECOLI,
    Constraints,
    OptimisationResult,
    build_table,
    optimise,
)
from gsynth_engine.duplex import DuplexView, construct_duplex
from gsynth_engine.align import PROTEIN_SCORING, Scoring, align, blosum62
from gsynth_engine.genbank import oligos_to_fasta, to_fasta, to_genbank
from gsynth_engine.ligation import ligation_series, plan_ligation
from gsynth_engine.primers import design_sequencing_primers
from gsynth_engine.verify import verify
from gsynth_engine.merzoug import AssemblyPlan, design_merzoug_assembly
from gsynth_engine.protocol import bench_protocol, order_sheet, order_sheet_csv
from gsynth_engine.sequence import SequenceError, gc_content
from gsynth_engine.ssd import SSDResult, design_small_sequence
from gsynth_engine.thermo import ANNEALING
from gsynth_engine import vectors as vector_catalogue
from rest_framework import status
from rest_framework.permissions import AllowAny
from rest_framework.throttling import ScopedRateThrottle
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.design.serializers import (
    CLEAVAGE_NAMES,
    CloneRequestSerializer,
    AlignRequestSerializer,
    LigationRequestSerializer,
    OptimiseRequestSerializer,
    PrimerRequestSerializer,
    VerifyRequestSerializer,
    resolve_vector,
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
        "warnings": result.warnings,
        # Empty means these two molecules really do join.
        "problems": result.problems,
        "is_clonable": result.is_clonable,
        "insert": _ssd_payload(ssd),
        "assembly": _assembly_payload(plan, result.name) if plan else None,
    }


def _spec_payload(spec) -> dict:
    """One catalogue entry, as the dropdown needs it."""
    return {
        "key": spec.key,
        "name": spec.name,
        "length": spec.length,
        "resistance": spec.resistance,
        "promoter": spec.promoter,
        "host": spec.host,
        "supplier": spec.supplier,
        "summary": spec.summary,
        "unique_sites": list(spec.unique_sites),
        "recommended_pairs": list(spec.recommended_pairs),
        "tags": [
            {"name": tag.name, "end": tag.end, "note": tag.note}
            for tag in spec.tags
        ],
        "notes": list(spec.notes),
        "reference": spec.reference,
        # Without a sequence the user has to import their own copy first.
        "has_sequence": spec.has_sequence,
        "supplies_translation_start": spec.supplies_translation_start,
        "tag_summary": spec.tag_summary,
    }


def _vector_payload(spec, sequence: str) -> dict:
    """Which vector was cut, and whether the sequence matches what it claims."""
    if spec is None:
        return {"recognised": False, "check": None, "spec": None}

    check = vector_catalogue.validate(sequence, spec)
    return {
        "recognised": True,
        "spec": _spec_payload(spec),
        "check": {
            "matches": check.matches,
            "length": check.length,
            "problems": check.problems,
            "notes": check.notes,
            "found_motifs": check.found_motifs,
            "missing_motifs": check.missing_motifs,
        },
    }


class VectorCatalogueView(APIView):
    """GET /api/design/vectors/ — the backbones G-Synth knows about.

    Public, like the enzyme table: it is reference data, and the cloning page
    needs it to build its dropdown before anything has been designed.
    """

    permission_classes = (AllowAny,)

    def get(self, request):
        return Response({
            "vectors": [_spec_payload(spec) for spec in vector_catalogue.CATALOGUE],
            "default": vector_catalogue.DEFAULT_VECTOR.key,
        })


class VectorSequenceView(APIView):
    """GET /api/design/vectors/<key>/ — a bundled sequence and its features."""

    permission_classes = (AllowAny,)

    def get(self, request, key: str):
        spec = vector_catalogue.get(key)
        if spec is None:
            return Response(
                {"detail": f"No vector called {key}."},
                status=status.HTTP_404_NOT_FOUND,
            )
        record = vector_catalogue.sequence_of(spec.key)
        if record is None:
            return Response(
                {"detail": f"{spec.name} has no bundled sequence. Import your "
                           f"own copy — G-Synth will check it against the entry."},
                status=status.HTTP_404_NOT_FOUND,
            )
        return Response({**record, "spec": _spec_payload(spec)})


def _optimisation_payload(result: OptimisationResult) -> dict:
    return {
        "sequence": result.sequence,
        "protein": result.protein,
        "length": result.length,
        "table": result.table,
        "cai_before": result.cai_before,
        "cai_after": result.cai_after,
        "gc_before": result.gc_before,
        "gc_after": result.gc_after,
        "sites_removed": result.sites_removed,
        "rare_codons_before": result.rare_codons_before,
        "rare_codons_after": result.rare_codons_after,
        "changed_codons": result.changed_codons,
        # Empty problems means the gene can be built and cut as asked.
        "problems": result.problems,
        "warnings": result.warnings,
        "is_clean": result.is_clean,
    }


class OptimiseView(APIView):
    """POST /api/design/optimise/ — rewrite a gene for the expression host.

    The protein is invariant; everything else is negotiable. Pass the cloning
    enzymes in `avoid_enzymes` so the result does not carry a site that would
    make the construct impossible to cut.
    """

    throttle_scope = "design"

    def post(self, request):
        serializer = OptimiseRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        table = ECOLI
        if data["reference_genes"]:
            try:
                table = build_table(
                    data["reference_genes"],
                    name="your reference set",
                    source=f"{len(data['reference_genes'])} genes you supplied",
                )
            except SequenceError as error:
                return _bad_request(error)

        constraints = Constraints(
            avoid_enzymes=tuple(data["avoid_enzymes"]),
            avoid_motifs=tuple(data["avoid_motifs"]),
            max_homopolymer=data["max_homopolymer"],
            gc_min=data["gc_min"],
            gc_max=data["gc_max"],
            gc_window=data["gc_window"],
            max_repeat=data["max_repeat"],
            avoid_rare=data["avoid_rare"],
        )
        try:
            result = optimise(
                data["sequence"],
                table=table,
                constraints=constraints,
                is_protein=data["is_protein"],
                keep_stop=data["keep_stop"],
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _optimisation_payload(result)
        payload["table_source"] = table.source
        return Response(payload)


class AlignView(APIView):
    """POST /api/design/align/ — compare two sequences.

    Separate from verification, which assumes the read is the construct and
    exploits that. This makes no such assumption: two genes from different
    strains, a design against what a supplier returned, a protein against
    its homologue.
    """

    throttle_scope = "design"

    def post(self, request):
        serializer = AlignRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        scoring = Scoring(
            match=data["match"],
            mismatch=data["mismatch"],
            gap_open=data["gap_open"],
            gap_extend=data["gap_extend"],
            matrix=blosum62() if data["is_protein"] else None,
        )
        try:
            result = align(
                data["first"], data["second"],
                mode=data["mode"], is_protein=data["is_protein"],
                scoring=scoring, try_reverse=data["try_reverse"],
            )
        except SequenceError as error:
            return _bad_request(error)

        return Response({
            "top": result.top,
            "marks": result.marks,
            "bottom": result.bottom,
            "rows": result.rows(60),
            "text": result.to_text(60),
            "score": result.score,
            "mode": result.mode,
            "length": result.length,
            "identity": result.identity,
            "similarity": result.similarity,
            "identities": result.identities,
            "similarities": result.similarities,
            "gaps": result.gaps,
            "start_a": result.start_a,
            "end_a": result.end_a,
            "start_b": result.start_b,
            "end_b": result.end_b,
            "reverse_complemented": result.reverse_complemented,
            "is_protein": result.is_protein,
            "warnings": result.warnings,
        })


class PrimerExportView(APIView):
    """POST /api/design/primers/export/?filetype=csv|fasta

    A primer set is ordered, not read on screen. CSV goes into a supplier's
    spreadsheet; FASTA into the ones that take an upload.
    """

    throttle_scope = "design"

    def post(self, request):
        serializer = PrimerRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        try:
            result = design_sequencing_primers(
                data["template"],
                target_start=data["target_start"], target_end=data["target_end"],
                circular=data["circular"], name=data["name"],
                tm_min=data["tm_min"], tm_max=data["tm_max"],
                margin=data["margin"], read_length=data["read_length"],
            )
        except SequenceError as error:
            return _bad_request(error)

        safe = (data["name"] or "seq").replace(" ", "_")
        if request.query_params.get("filetype") == "fasta":
            return _attachment(
                oligos_to_fasta(result.as_rows),
                f"{safe}_primers.fasta", "text/plain; charset=utf-8",
            )

        import csv
        import io

        buffer = io.StringIO()
        rows = result.as_rows
        if rows:
            writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        return _attachment(buffer.getvalue(), f"{safe}_primers.csv", "text/csv")


class LigationView(APIView):
    """POST /api/design/ligation/ — masses to pipette, at several ratios.

    Nobody runs one ligation, so the response is the whole series.
    """

    def post(self, request):
        serializer = LigationRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        common = {
            "vector_length": data["vector_length"],
            "insert_length": data["insert_length"],
            "vector_ng": data["vector_ng"],
            "ends": data["ends"],
        }
        try:
            plans = (
                ligation_series(**common, ratios=tuple(data["ratios"]))
                if data["ratios"]
                else [plan_ligation(**common, total_volume_uL=data["total_volume_uL"])]
            )
        except SequenceError as error:
            return _bad_request(error)

        return Response({
            "reactions": [
                {
                    "ratio": plan.ratio,
                    "vector_ng": plan.vector_ng,
                    "insert_ng": plan.insert_ng,
                    "vector_fmol": plan.vector_fmol,
                    "insert_fmol": plan.insert_fmol,
                    "total_ng": plan.total_ng,
                    "rows": plan.as_rows(),
                    "warnings": plan.warnings,
                }
                for plan in plans
            ],
            "vector_length": data["vector_length"],
            "insert_length": data["insert_length"],
            "ends": data["ends"],
            "total_volume_uL": data["total_volume_uL"],
        })


class SequencingPrimerView(APIView):
    """POST /api/design/primers/ — primers that read across a region."""

    throttle_scope = "design"

    def post(self, request):
        serializer = PrimerRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        try:
            result = design_sequencing_primers(
                data["template"],
                target_start=data["target_start"],
                target_end=data["target_end"],
                circular=data["circular"],
                name=data["name"],
                tm_min=data["tm_min"],
                tm_max=data["tm_max"],
                margin=data["margin"],
                read_length=data["read_length"],
            )
        except SequenceError as error:
            return _bad_request(error)

        return Response({
            "primers": [
                {
                    "name": p.name,
                    "sequence": p.sequence,
                    "length": p.length,
                    "start": p.start,
                    "direction": p.direction,
                    "tm": p.tm,
                    "gc": p.gc,
                    "reads_from": p.reads_from,
                    "reads_to": p.reads_to,
                }
                for p in result.primers
            ],
            "rows": result.as_rows,
            "target_start": result.target_start,
            "target_end": result.target_end,
            "gaps": result.gaps,
            "covers_target": result.covers_target,
            "warnings": result.warnings,
        })


class VerifyView(APIView):
    """POST /api/design/verify/ — do the reads say you built the design?"""

    throttle_scope = "design"

    def post(self, request):
        serializer = VerifyRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        region = None
        if data.get("region_start") is not None and data.get("region_end") is not None:
            region = (data["region_start"], data["region_end"])

        try:
            report = verify(
                data["design"], data["reads"],
                circular=data["circular"], trim=data["trim"],
                coding_start=data.get("coding_start"),
                coding_end=data.get("coding_end"),
                region=region,
            )
        except SequenceError as error:
            return _bad_request(error)

        return Response({
            "design_length": report.design_length,
            "coverage": report.coverage,
            "gaps": report.gaps,
            "fully_covered": report.fully_covered,
            # Empty differences with at least one read means it is the design.
            "is_verified": report.is_verified,
            "differences": [
                {
                    "kind": d.kind,
                    "position": d.position,
                    "expected": d.expected,
                    "found": d.found,
                    "residue": d.residue,
                    "from_residue": d.from_residue,
                    "to_residue": d.to_residue,
                    "silent": d.silent,
                    "description": d.description,
                }
                for d in report.differences
            ],
            "reads": [
                {
                    "name": r.name,
                    "length": r.length,
                    "start": r.start,
                    "end": r.end,
                    "covered": r.covered,
                    "reverse_complemented": r.reverse_complemented,
                    "identity": r.identity,
                    "matched": r.matched,
                    "difference_count": len(r.differences),
                    "is_clean": r.is_clean,
                    "warnings": r.warnings,
                }
                for r in report.reads
            ],
            "warnings": report.warnings,
        })


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

    throttle_scope = "design"

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

    throttle_scope = "design"

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

    throttle_scope = "design"

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

            vector_sequence, vector_name, annotations, spec = resolve_vector(data)
            result = clone(
                vector_sequence,
                insert_forward,
                insert_reverse=insert_reverse,
                left_enzyme=data["left_enzyme"],
                right_enzyme=data["right_enzyme"],
                circular=data["vector_is_circular"],
                name=name,
                vector_annotations=annotations,
                vector_spec=spec,
                orf_start=ssd.orf_start,
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _clone_payload(result, ssd, plan)
        payload["vector_name"] = vector_name
        payload["vector"] = _vector_payload(spec, vector_sequence)

        if data["save_as_project"] and result.is_clonable:
            project = Project.objects.create(
                user=request.user,
                name=name,
                module="cloning",
                sequence=result.plasmid,
                notes=f"{result.length} bp in {vector_name} · "
                      f"{result.left_enzyme}/{result.right_enzyme}",
                data=payload,
            )
            payload["project_id"] = project.id
        return Response(payload)


def _attachment(text: str, filename: str, content_type: str) -> HttpResponse:
    response = HttpResponse(text, content_type=content_type)
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response


def _today() -> str:
    """GenBank's date field. Taken here rather than in the engine, so the
    engine's own output stays byte-identical between runs."""
    from django.utils import timezone

    return timezone.now().strftime("%d-%b-%Y").upper()


class CloneExportView(APIView):
    """POST /api/design/clone/export/?filetype=genbank|fasta — as a file.

    GenBank by default, because that is what carries the features across to
    SnapGene, Benchling or ApE. A design that only exists inside G-Synth is
    not finished.
    """

    throttle_scope = "design"

    def post(self, request):
        serializer = CloneRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        name = data["name"]

        try:
            vector_sequence, vector_name, annotations, spec = resolve_vector(data)
            if data["fragment"]:
                plan = design_merzoug_assembly(
                    data["sequence"], **serializer.engine_kwargs
                )
                ssd, forward, reverse = plan.ssd, plan.construct_forward, plan.construct_reverse
            else:
                kwargs = {
                    k: v for k, v in serializer.engine_kwargs.items()
                    if k not in ("target_oligo_length", "overhang_length")
                }
                ssd = design_small_sequence(data["sequence"], **kwargs)
                forward, reverse = ssd.forward, ssd.reverse

            result = clone(
                vector_sequence, forward, insert_reverse=reverse,
                left_enzyme=data["left_enzyme"], right_enzyme=data["right_enzyme"],
                circular=data["vector_is_circular"], name=name,
                vector_annotations=annotations, vector_spec=spec,
                orf_start=ssd.orf_start,
            )
        except SequenceError as error:
            return _bad_request(error)

        safe = (name or "construct").replace(" ", "_")
        if request.query_params.get("filetype") == "fasta":
            return _attachment(
                to_fasta(
                    result.plasmid, name=safe,
                    description=f"{result.length} bp in {vector_name}",
                ),
                f"{safe}.fasta", "text/plain; charset=utf-8",
            )

        payload = _clone_payload(result, ssd, None)
        return _attachment(
            to_genbank(
                result.plasmid,
                name=safe,
                description=f"{name} cloned into {vector_name} "
                            f"({result.left_enzyme}/{result.right_enzyme})",
                features=payload["annotations"],
                circular=True,
                date=_today(),
            ),
            f"{safe}.gb", "chemical/seq-na-genbank",
        )


class ConstructExportView(APIView):
    """POST /api/design/assembly/export/ — the construct, or its oligos.

    `filetype=oligos` gives one FASTA entry per oligo, which is what a supplier
    accepts as an upload. Retyping thirty oligo names into a web form is
    where transcription errors come from.
    """

    throttle_scope = "design"

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

        safe = (name or "construct").replace(" ", "_")
        wanted = request.query_params.get("filetype", "genbank")

        if wanted == "oligos":
            return _attachment(
                oligos_to_fasta(
                    [o.as_row for o in order_sheet(plan, construct_name=name)]
                ),
                f"{safe}_oligos.fasta", "text/plain; charset=utf-8",
            )
        if wanted == "fasta":
            return _attachment(
                to_fasta(plan.construct_forward, name=safe),
                f"{safe}.fasta", "text/plain; charset=utf-8",
            )

        features = [
            {"name": s.name, "type": "misc_feature", "start": s.start,
             "end": s.end, "direction": 1}
            for s in plan.ssd.segments
        ] + [
            {"name": f.name, "type": "misc_feature", "start": f.top_start,
             "end": f.top_end, "direction": 1}
            for f in plan.fragments
        ]
        return _attachment(
            to_genbank(
                plan.construct_forward, name=safe,
                description=f"{name}: {plan.fragment_count} fragments, "
                            f"{plan.oligo_count} oligos",
                features=features, circular=False, date=_today(),
            ),
            f"{safe}.gb", "chemical/seq-na-genbank",
        )


class OrderSheetView(APIView):
    """POST /api/design/assembly/order-sheet/ — the oligo list as CSV."""

    throttle_scope = "design"

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

    throttle_scope = "design"

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
