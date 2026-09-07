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
    FeatureDetectionSerializer,
    HybridizationRequestSerializer,
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
from apps.sequences.sbol import to_sbol3
from gsynth_engine import vectors as vector_catalogue
from gsynth_engine.align import Scoring, align, blosum62
from gsynth_engine.annotations import detect_common_features
from gsynth_engine.chromatogram import read_trace, summarise
from gsynth_engine.cloning import (
    CloningResult,
    clone,
    find_sites,
    open_reading_frames,
)
from gsynth_engine.codon import (
    CODON_DATA_SHA256,
    CODON_DATA_URL,
    CODON_DATA_VERSION,
    DEFAULT_HOST,
    TABLES,
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
from gsynth_engine.esd import ESDResult, design_extended_sequence
from gsynth_engine.gel import ladder_payload, recommended_ladder, restriction_digest_sizes
from gsynth_engine.genbank import oligos_to_fasta, to_fasta, to_genbank
from gsynth_engine.hybridization import hybridize
from gsynth_engine.ligation import ligation_series, plan_ligation
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
from gsynth_engine.thermo import ANNEALING, BufferConditions
from gsynth_engine.verify import ConsensusReport, assemble_consensus, verify


def _bad_request(error: SequenceError) -> Response:
    """Engine errors are already written for the user — pass them through."""
    return Response({"detail": str(error)}, status=status.HTTP_400_BAD_REQUEST)


def _gel_simulation(title: str, lanes: list[dict], fragment_sizes: list[int]) -> dict:
    """A plotting-ready prediction, explicitly not an experimental image."""
    return {
        "title": title,
        "prediction_only": True,
        "notice": (
            "In-silico size prediction. Migration, intensity, topology effects, "
            "partial digestion and background require experimental confirmation."
        ),
        "recommended_ladder": recommended_ladder(fragment_sizes),
        "ladders": ladder_payload(),
        "lanes": lanes,
    }


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


def _assembly_payload(plan: ESDResult, construct_name: str) -> dict:
    insert = next(segment for segment in plan.ssd.segments if segment.name == "insert")
    return {
        "construct_forward": plan.construct_forward,
        "construct_reverse": plan.construct_reverse,
        "construct_length": plan.construct_length,
        "construct_gc": round(gc_content(plan.construct_forward), 1),
        # Zero-based, end-exclusive coordinates let the saved linear design
        # flow directly into sequencing-primer and read verification tools.
        # There is deliberately no backbone_length here: this is an assembly
        # cassette, not yet a cloned plasmid, so a ligation mass calculation
        # would otherwise mistake the cassette for the vector backbone.
        "insert_start": insert.start,
        "insert_end": insert.end,
        "topology": "linear",
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
    """Every catalogued restriction occurrence as a filterable map feature.

    The client defaults to single-cutters plus the cloning pair because every
    four-base cutter at once is visual noise.  Returning the complete set is
    still essential: an explicit “all sites” choice must actually show all
    sites, rather than pretending a multi-cutter does not exist.
    """
    out: list[dict] = []
    for enzyme in sorted(ALL_ENZYMES):
        sites = find_sites(plasmid, enzyme, circular=True)
        if not sites:
            continue
        used = enzyme in highlight
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
    frame = result.reading_frame

    def check(label: str, status: str, detail: str) -> dict:
        return {
            "check": label,
            "status": status,
            "passed": status == "pass",
            "detail": detail,
        }

    return [
        check(
            "Overhangs are compatible",
            "pass" if ends_ok else "block",
            " ".join(p for p in result.problems if "does not match" in p)
            or f"{result.left_enzyme} and {result.right_enzyme} ends match the vector.",
        ),
        check(
            "Both strands pair everywhere",
            "pass" if not duplex_mismatches and not duplex_problems else "block",
            (
                f"{len(duplex_mismatches)} positions do not pair."
                if duplex_mismatches else " ".join(duplex_problems)
                if duplex_problems else "No mismatch in the insert duplex."
            ),
        ),
        check(
            "Each enzyme cuts the vector once",
            "pass",
            "Checked when the vector was cut; a second site is refused.",
        ),
        check(
            "Orientation is forced",
            "pass" if result.left_enzyme != result.right_enzyme else "block",
            (
                "The insert reads on the minus strand of the vector's own "
                "numbering, as every pET cassette does."
                if result.reversed_insert
                else "The insert reads with the vector's numbering."
            ),
        ),
        check(
            "Sites are regenerated at both seams",
            "pass" if junctions_ok else "review",
            (
                "The insert can be cut back out — which is how the clone gets "
                "verified on a gel."
                if junctions_ok else
                "At least one site is lost, so the insert cannot be excised."
            ),
        ),
        check(
            "Expression reading frame",
            "pass" if frame.status == "not_applicable" else frame.status,
            frame.summary,
        ),
    ]


def _clone_payload(result: CloningResult, ssd: SSDResult | None, plan, insert_annotations=None) -> dict:
    """The recombinant plasmid, shaped for a map viewer.

    The insert is sent as one more annotation so the client can draw the
    whole plasmid from a single list, rather than special-casing it.
    """
    annotations = list(result.annotations)
    cassette_annotation = {
        "name": result.name,
        "type": "CDS" if result.protein else "misc_feature",
        "start": result.insert_start,
        "end": result.insert_end,
        "direction": 1,
        "color": "#0E6E77",
    }
    # Duplex handoffs have no SSD object. Use the same origin as the engine's
    # protein so cohesive-end bases cannot shift the displayed codon frame.
    if result.translation_start is not None and result.protein:
        cassette_annotation["translation_start"] = result.translation_start
        cassette_annotation["translation_end"] = result.insert_end
    annotations.append(cassette_annotation)
    for annotation in insert_annotations or []:
        moved = {**annotation, "start": result.insert_start + annotation["start"],
                 "end": result.insert_start + annotation["end"]}
        for key in ("translation_start", "translation_end"):
            if annotation.get(key) is not None:
                moved[key] = result.insert_start + annotation[key]
        annotations.append(moved)

    # Preserve the cassette's internal biological meaning on the project map.
    # These coordinates come from the same SSD result that built the oligo;
    # they are not re-detected from motifs, so a repeated His or linker motif
    # cannot make an annotation jump to the wrong occurrence.
    segment_colours = {
        "overhang": "#C97634",
        "start codon": "#9E3D3D",
        "linker": "#78889B",
        "6×his tag": "#0E6E77",
        "site": "#6A4C93",
        "insert": "#3F7A52",
    }
    if ssd is not None:
        for segment in ssd.segments:
            key = next(
                (name for name in segment_colours if name in segment.name.lower()),
                "linker",
            )
            label = (
                f"{result.name} target"
                if segment.name.lower() == "insert"
                else segment.name
            )
            annotations.append({
                "name": label,
                "type": "mat_peptide" if segment.name.lower() == "insert" else "misc_feature",
                "start": result.insert_start + segment.start,
                "end": result.insert_start + segment.end,
                "direction": 1,
                "color": segment_colours[key],
            })

    duplex_mismatches = construct_duplex(plan).mismatches() if plan else []
    try:
        digest_sizes = restriction_digest_sizes(
            result.plasmid,
            [result.left_enzyme, result.right_enzyme],
            circular=True,
        )
        digest_lane = {
            "name": f"{result.left_enzyme} + {result.right_enzyme}",
            "description": "Complete diagnostic double digest",
            "bands": [{"size_bp": size, "label": f"{size:,} bp"} for size in digest_sizes],
        }
        gel = _gel_simulation(
            "Predicted diagnostic digest",
            [digest_lane],
            digest_sizes,
        )
    except SequenceError as error:
        gel = _gel_simulation("Predicted diagnostic digest", [], [])
        gel["notice"] = f"Diagnostic digest unavailable: {error}"

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
        "reading_frame": {
            "status": result.reading_frame.status,
            "confirmed": result.reading_frame.confirmed,
            "summary": result.reading_frame.summary,
            "translation_start": result.reading_frame.translation_start,
            "start_codon": result.reading_frame.start_codon,
            "start_source": result.reading_frame.start_source,
            "rbs_name": result.reading_frame.rbs_name,
            "rbs_start": result.reading_frame.rbs_start,
            "rbs_end": result.reading_frame.rbs_end,
            "rbs_spacing_nt": result.reading_frame.rbs_spacing_nt,
            "rbs_source": result.reading_frame.rbs_source,
            "promoter_name": result.reading_frame.promoter_name,
            "promoter_start": result.reading_frame.promoter_start,
            "promoter_source": result.reading_frame.promoter_source,
            "stop_codon": result.reading_frame.stop_codon,
            "stop_position": result.reading_frame.stop_position,
            "stop_context": result.reading_frame.stop_context,
            "left_junction_offset": result.reading_frame.left_junction_offset,
            "right_junction_phase": result.reading_frame.right_junction_phase,
            "protein_length": result.reading_frame.protein_length,
            "checks": [
                {
                    "code": frame_check.code,
                    "label": frame_check.label,
                    "status": frame_check.status,
                    "detail": frame_check.detail,
                }
                for frame_check in result.reading_frame.checks
            ],
        },
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
        "gel": gel,
        "validation": _validation(
            result,
            duplex_mismatches,
        ),
        "preflight": cloning_preflight(
            result, duplex_mismatches=duplex_mismatches,
        ).to_dict(),
        "warnings": result.warnings,
        # Empty when all cloning checks pass.
        "problems": result.problems,
        "is_clonable": result.is_clonable,
        # Pre-digested inserts have no generated SSD payload.
        "insert": _ssd_payload(ssd) if ssd is not None else None,
        "assembly": _assembly_payload(plan, result.name) if plan else None,
    }


def _reviewed_product_annotations(data: dict, sequence_length: int) -> list[dict] | None:
    """Validate optional user-reviewed features against the product molecule.

    The request serializer checks each feature's shape. Product length only
    exists after cloning, so coordinate bounds are enforced here before the
    list can enter a saved project or GenBank export. Circular origin-crossing
    features may extend once beyond ``sequence_length``.
    """
    annotations = data.get("product_annotations")
    if annotations is None:
        return None
    reviewed = [dict(annotation) for annotation in annotations]
    for annotation in reviewed:
        start = annotation["start"]
        end = annotation["end"]
        if start >= sequence_length or end > sequence_length * 2 or end - start > sequence_length:
            raise SequenceError(
                f"Feature ‘{annotation['name']}’ lies outside the {sequence_length:,} bp product."
            )
        if end > sequence_length and end - sequence_length >= start:
            raise SequenceError(
                f"Feature ‘{annotation['name']}’ does not cross the circular origin correctly."
            )
    return reviewed


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
        "expression_capable": spec.expression_capable,
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
        "input_protein": result.input_protein,
        "protein_context": result.protein_context,
        "initiator_methionine_added": result.initiator_methionine_added,
        "recommended_design_is_coding": result.recommended_design_is_coding,
    }


class CodonHostCatalogueView(APIView):
    """GET /api/design/codon-hosts/ — reproducible bundled host profiles."""

    permission_classes = (AllowAny,)

    def get(self, request):
        return Response({
            "default": DEFAULT_HOST,
            "dataset": {
                "name": "FDA HIVE-CUTs / CoCoPUTs",
                "release": CODON_DATA_VERSION,
                "url": CODON_DATA_URL,
                "sha256": CODON_DATA_SHA256,
            },
            "hosts": [
                {
                    "key": key,
                    "name": table.name,
                    "source": table.source,
                    "category": table.category,
                    "taxon_id": table.taxon_id,
                    "dataset": table.dataset,
                    "dataset_release": table.dataset_release,
                    "data_scope": table.data_scope,
                    "coding_sequences": table.coding_sequences,
                    "codon_count": table.codon_count,
                    "gc_percent": table.gc_percent,
                    "source_url": table.source_url,
                    "metric_label": "Profile-relative CAI",
                }
                for key, table in TABLES.items()
            ],
        })


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

        table = TABLES[data["host"]]
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
                protein_context=data["protein_context"],
                keep_stop=data["keep_stop"],
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _optimisation_payload(result)
        payload["host"] = data["host"] if not data["reference_genes"] else "custom"
        payload["table_source"] = table.source
        payload["metric_label"] = (
            "CAI (custom reference set)"
            if data["reference_genes"]
            else "Profile-relative CAI"
        )
        payload["expression_yield_predicted"] = False
        payload["preflight"] = optimisation_preflight(result).to_dict()
        payload["provenance"] = _provenance(
            "codon_optimisation", data, result.sequence,
        )
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


class HybridizationView(APIView):
    """POST /api/design/hybridize/ — physical antiparallel strand pairing."""

    throttle_scope = "design"

    def post(self, request):
        serializer = HybridizationRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        conditions = BufferConditions(
            name="hybridization analysis",
            oligo_nM=data["oligo_nM"],
            na_mM=data["na_mM"],
            mg_mM=data["mg_mM"],
            dntp_mM=data["dntp_mM"],
        )
        try:
            result = hybridize(
                data["first"],
                data["second"],
                conditions=conditions,
                analysis_temperature_c=data["analysis_temperature_c"],
            )
        except SequenceError as error:
            return _bad_request(error)

        return Response({
            "first": result.first,
            "second": result.second,
            "top": result.top,
            "marks": result.marks,
            "bottom": result.bottom,
            "rows": result.rows(60),
            "width": result.width,
            "offset": result.offset,
            "overlap_start": result.overlap_start,
            "overlap_end": result.overlap_end,
            "overlap_length": result.overlap_length,
            "paired_bases": result.paired_bases,
            "paired_percent": result.paired_percent,
            "mismatches": result.mismatches,
            "longest_perfect_run": result.longest_perfect_run,
            "complementarity": result.complementarity,
            "predicted_state": result.predicted_state,
            "overhangs": [item.payload() for item in result.overhangs],
            "left_end": result.left_end,
            "right_end": result.right_end,
            "alternative_placements": result.alternative_placements,
            "tm_c": result.tm_c,
            "tm_margin_c": result.tm_margin_c,
            "delta_h_kcal_mol": result.delta_h_kcal_mol,
            "delta_s_cal_mol_k": result.delta_s_cal_mol_k,
            "analysis_temperature_c": result.analysis_temperature_c,
            "conditions": {
                "name": result.conditions.name,
                "oligo_nM": result.conditions.oligo_nM,
                "na_mM": result.conditions.na_mM,
                "mg_mM": result.conditions.mg_mM,
                "dntp_mM": result.conditions.dntp_mM,
                "summary": result.conditions.summary,
            },
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
    """Return a practical series of insert-to-vector ligation ratios."""

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


def _difference_payload(d) -> dict:
    return {
        "kind": d.kind,
        "position": d.position,
        "expected": d.expected,
        "found": d.found,
        "residue": d.residue,
        "from_residue": d.from_residue,
        "to_residue": d.to_residue,
        "silent": d.silent,
        "description": d.description,
        # None when the read came as letters. False means the trace does not
        # support it — noise, not a mutation.
        "quality": d.quality,
        "confident": d.confident,
        "read_index": d.read_index,
    }


def _read_payload(r) -> dict:
    return {
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
        "mean_quality": r.mean_quality,
        "trimmed_start": r.trimmed_start,
        "trimmed_end": r.trimmed_end,
    }


def _consensus_payload(report: ConsensusReport) -> dict:
    return {
        "sequence": report.sequence,
        "coverage": report.coverage,
        "identity": report.identity,
        "fully_covered": report.fully_covered,
        "bidirectional_overlap": report.bidirectional_overlap,
        "bidirectional_agreement": report.bidirectional_agreement,
        "gaps": report.gaps,
        "difference_count": len(report.differences),
        "warnings": report.warnings,
    }


class TraceVerifyView(APIView):
    """POST /api/design/verify/traces/ — the reads, with their peaks.

    The same comparison as `/verify/`, except the reads arrive as ABIF or SCF
    chromatogram files.
    That buys two things the letters cannot give: the ends are trimmed by
    quality rather than by a fixed count, and every difference is returned
    with the confidence of the base that produced it — plus the slice of
    trace around it, so the drawing can show the peak that was called.
    """

    throttle_scope = "design"
    parser_classes = [MultiPartParser, FormParser]

    def post(self, request):
        serializer = TraceUploadSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        traces, summaries = {}, []
        for upload in data["traces"]:
            name = upload.name or f"trace {len(traces) + 1}"
            try:
                trace = read_trace(upload.read(), name=name)
            except SequenceError as error:
                return _bad_request(error)
            traces[name] = trace
            s = summarise(trace, data["trim_quality"])
            summaries.append({
                "name": s.name, "length": s.length,
                "mean_quality": s.mean_quality,
                "trim_start": s.trim_start, "trim_stop": s.trim_stop,
                "trimmed_length": s.trimmed_length,
                "high_quality_bases": s.high_quality_bases,
                "sample_count": s.sample_count, "usable": s.usable,
            })

        region = None
        if data.get("region_start") is not None and data.get("region_end") is not None:
            region = (data["region_start"], data["region_end"])

        try:
            report = verify(
                data["design"], {n: t.sequence for n, t in traces.items()},
                circular=data["circular"],
                coding_start=data.get("coding_start"),
                coding_end=data.get("coding_end"),
                region=region, traces=traces,
                trim_quality=data["trim_quality"],
            )
            consensus = assemble_consensus(
                data["design"], traces,
                circular=data["circular"], region=region,
                trim_quality=data["trim_quality"],
            )
            raw_consensus = assemble_consensus(
                data["design"], traces,
                circular=data["circular"], region=region,
                trim_quality=0,
            )
        except SequenceError as error:
            return _bad_request(error)

        # The peaks around each difference, so it can be looked at rather
        # than taken on trust. Only these windows travel, never whole traces.
        windows = []
        tracks = []
        for read in report.reads:
            trace = traces.get(read.name)
            if trace is None:
                continue
            track = trace.alignment_track(
                read.trimmed_start,
                trace.length - read.trimmed_end,
                reverse=read.reverse_complemented,
            )
            tracks.append({
                "read": read.name,
                "reference_start": read.start,
                "reference_end": read.end,
                "reverse_complemented": read.reverse_complemented,
                **track,
            })
            for d in read.differences:
                if d.read_index is None:
                    continue
                windows.append({
                    "read": read.name,
                    "position": d.position,
                    **trace.window(d.read_index),
                })

        preflight = verification_preflight(report)
        return Response({
            "design_length": report.design_length,
            "region_start": region[0] if region else 0,
            "region_end": region[1] if region else report.design_length,
            "coverage": report.coverage,
            "quality_cutoff": data["trim_quality"],
            "consensus": _consensus_payload(consensus),
            "raw_consensus": _consensus_payload(raw_consensus),
            "gaps": report.gaps,
            "fully_covered": report.fully_covered,
            "is_verified": report.is_verified,
            "verification_state": verification_state(report),
            "preflight": preflight.to_dict(),
            "provenance": _provenance("sequence_verification", data, data["design"]),
            "differences": [_difference_payload(d) for d in report.differences],
            "reads": [_read_payload(r) for r in report.reads],
            "traces": summaries,
            "trace_tracks": tracks,
            "trace_windows": windows,
            "warnings": report.warnings,
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

        preflight = verification_preflight(report)
        return Response({
            "design_length": report.design_length,
            "region_start": region[0] if region else 0,
            "region_end": region[1] if region else report.design_length,
            "coverage": report.coverage,
            "gaps": report.gaps,
            "fully_covered": report.fully_covered,
            # Verification requires agreement over the complete requested region.
            "is_verified": report.is_verified,
            "verification_state": verification_state(report),
            "preflight": preflight.to_dict(),
            "provenance": _provenance("sequence_verification", data, data["design"]),
            "differences": [_difference_payload(d) for d in report.differences],
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
        for name in sorted(ALL_ENZYMES):
            sequence, kind = overhang(name)
            aliases = list(ALL_ENZYMES[name].get("aliases", ()))
            enzymes.append({
                "name": name,
                "aliases": aliases,
                "recognition": ALL_ENZYMES[name]["recognition"],
                # Preferred cloning enzymes are offered first.
                "common": name in RESTRICTION_ENZYMES,
                "overhang": sequence,
                "overhang_type": kind,
                # This is derived from the retained top-strand remainder.
                # An ATG elsewhere in the recognition site may be cut away
                # or followed by frame-shifting bases.
                "supplies_start_codon": supplies_start_codon(name),
            })
        return Response({
            "enzymes": enzymes,
            "canonical_geometries": len(enzymes),
            "selectable_names": sum(1 + len(enzyme["aliases"]) for enzyme in enzymes),
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
        payload["provenance"] = _provenance(
            "ssd", serializer.validated_data, result.forward,
        )

        if serializer.validated_data["save_as_project"]:
            project = Project.objects.create(
                user=request.user,
                name=serializer.validated_data["name"],
                module="ssd",
                sequence=result.forward,
                notes=f"{result.left_enzyme}/{result.right_enzyme}"
                      + (f" · {result.cleavage_site}" if result.cleavage_site else ""),
                data=payload,
                provenance=payload["provenance"],
            )
            payload["project_id"] = project.id
        return Response(payload)


class ExtendedSequenceDesignView(APIView):
    """POST /api/design/assembly/ — one insert, an ordered set of oligo pairs."""

    throttle_scope = "design"

    def post(self, request):
        serializer = SaveableAssemblyRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        name = serializer.validated_data["name"]
        try:
            plan = design_extended_sequence(
                serializer.validated_data["sequence"], **serializer.engine_kwargs
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _assembly_payload(plan, name)
        payload["provenance"] = _provenance(
            "extended_sequence_design", serializer.validated_data, plan.construct_forward,
        )

        if serializer.validated_data["save_as_project"]:
            project = Project.objects.create(
                user=request.user,
                name=name,
                module="extended_sequence_design",
                sequence=plan.construct_forward,
                notes=f"{plan.fragment_count} fragments · "
                      f"{plan.oligo_count} oligos · "
                      f"{plan.overhang_length} nt overhangs",
                data=payload,
                provenance=payload["provenance"],
            )
            payload["project_id"] = project.id
        return Response(payload)


def _run_clone(data: dict, engine_kwargs: dict):
    """Build exactly one cloning result for preview, export and worksheet."""
    if data.get("pre_digested"):
        plan = None
        ssd = None
        insert_forward = data["sequence"]
        insert_reverse = data["insert_reverse"]
    elif data["fragment"]:
        plan = design_extended_sequence(data["sequence"], **engine_kwargs)
        ssd = plan.ssd
        insert_forward = plan.construct_forward
        insert_reverse = plan.construct_reverse
    else:
        plan = None
        ssd = design_small_sequence(
            data["sequence"],
            **{key: value for key, value in engine_kwargs.items()
               if key not in ("target_oligo_length", "overhang_length")},
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
        name=data["name"],
        vector_annotations=annotations,
        vector_spec=spec,
        orf_start=(
            data.get("orf_start")
            if data.get("pre_digested")
            else ssd.orf_start if ssd is not None else None
        ),
        auto_detect_frame=bool(data.get("pre_digested")),
    )
    return result, ssd, plan, vector_sequence, vector_name, spec


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
            result, ssd, plan, vector_sequence, vector_name, spec = _run_clone(
                data, serializer.engine_kwargs,
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = _clone_payload(result, ssd, plan, data.get("insert_annotations"))
        try:
            reviewed_annotations = _reviewed_product_annotations(data, result.length)
        except SequenceError as error:
            return _bad_request(error)
        if reviewed_annotations is not None:
            payload["annotations"] = reviewed_annotations
        payload["vector_name"] = vector_name
        payload["vector"] = _vector_payload(spec, vector_sequence)
        payload["provenance"] = _provenance(
            "cloning", data, result.plasmid, vector=vector_sequence,
        )

        if data["save_as_project"] and result.is_clonable:
            project = Project.objects.create(
                user=request.user,
                name=name,
                module="cloning",
                sequence=result.plasmid,
                notes=f"{result.length} bp in {vector_name} · "
                      f"{result.left_enzyme}/{result.right_enzyme}",
                data=payload,
                provenance=payload["provenance"],
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
    """POST /api/design/clone/export/?filetype=genbank|fasta|sbol3."""

    throttle_scope = "design"

    def post(self, request):
        serializer = CloneRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        name = data["name"]

        try:
            result, ssd, plan, vector_sequence, vector_name, _spec = _run_clone(
                data, serializer.engine_kwargs,
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

        payload = _clone_payload(result, ssd, plan, data.get("insert_annotations"))
        try:
            reviewed_annotations = _reviewed_product_annotations(data, result.length)
        except SequenceError as error:
            return _bad_request(error)
        if reviewed_annotations is not None:
            payload["annotations"] = reviewed_annotations
        if request.query_params.get("filetype") == "sbol3":
            return _attachment(
                to_sbol3(
                    result.plasmid,
                    name=name,
                    description=f"{name} cloned into {vector_name} "
                                f"({result.left_enzyme}/{result.right_enzyme})",
                    features=payload["annotations"],
                    circular=True,
                ),
                f"{safe}.sbol.json", "application/ld+json; charset=utf-8",
            )
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


class CloneWorksheetView(APIView):
    """POST /api/design/clone/worksheet/ — printable release-to-bench record."""

    throttle_scope = "design"

    def post(self, request):
        serializer = CloneRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        try:
            result, _ssd, _plan, vector_sequence, vector_name, _spec = _run_clone(
                data, serializer.engine_kwargs,
            )
            if not result.is_clonable:
                raise SequenceError(
                    "A bench worksheet cannot be released for a blocked cloning design."
                )
            primers = design_sequencing_primers(
                result.plasmid,
                target_start=result.insert_start,
                target_end=result.insert_end,
                circular=True,
                name=(data["name"] or "clone").replace(" ", "_")[:20],
            )
            reactions = ligation_series(
                vector_length=result.backbone_length,
                insert_length=result.insert_length,
                ends=result.junctions[0].kind if result.junctions else "5'",
            )
        except SequenceError as error:
            return _bad_request(error)

        preflight = cloning_preflight(result)
        provenance = _provenance(
            "cloning", data, result.plasmid, vector=vector_sequence,
        )
        text = cloning_worksheet(
            result,
            vector_name=vector_name,
            primer_set=primers,
            ligation_plans=reactions,
            preflight=preflight,
            provenance=provenance,
        )
        safe = (data["name"] or "construct").replace(" ", "_")
        return _attachment(
            text, f"{safe}_bench_worksheet.txt", "text/plain; charset=utf-8",
        )


class ConstructExportView(APIView):
    """POST /api/design/assembly/export/ — construct and oligo sequences.

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
            plan = design_extended_sequence(
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
        if wanted == "all-sequences":
            sequence_set = "".join([
                to_fasta(
                    plan.construct_forward,
                    name=f"{safe}_assembled_forward",
                    description="assembled construct, forward strand, 5-prime to 3-prime",
                ),
                to_fasta(
                    plan.construct_reverse,
                    name=f"{safe}_assembled_reverse",
                    description="assembled construct, reverse strand, 5-prime to 3-prime",
                ),
                oligos_to_fasta(
                    [o.as_row for o in order_sheet(plan, construct_name=name)]
                ),
            ])
            return _attachment(
                sequence_set,
                f"{safe}_all_sequences.fasta", "text/plain; charset=utf-8",
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
        if wanted == "sbol3":
            return _attachment(
                to_sbol3(
                    plan.construct_forward,
                    name=name,
                    description=f"{name}: {plan.fragment_count} fragments, "
                                f"{plan.oligo_count} oligos",
                    features=features,
                    circular=False,
                ),
                f"{safe}.sbol.json", "application/ld+json; charset=utf-8",
            )
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
            plan = design_extended_sequence(
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
            plan = design_extended_sequence(
                serializer.validated_data["sequence"], **serializer.engine_kwargs
            )
        except SequenceError as error:
            return _bad_request(error)

        text = bench_protocol(plan, construct_name=name)
        response = HttpResponse(text, content_type="text/plain; charset=utf-8")
        safe_name = name.replace(" ", "_") or "construct"
        response["Content-Disposition"] = f'attachment; filename="{safe_name}_protocol.txt"'
        return response


def _primer_payload(primer) -> dict:
    """One primer, with both Tm figures kept distinct.

    `tm` is the annealing portion and is what the annealing temperature was
    derived from; `tm_full` is the whole oligo and only applies once the tail
    has been copied. Collapsing them into one number is what makes a tailed
    primer look like it should anneal ten degrees hotter than it does.
    """
    return {
        "name": primer.name,
        "sequence": primer.sequence,
        "tail": primer.tail,
        "anneals": primer.anneals,
        "direction": primer.direction,
        "start": primer.start,
        "end": primer.end,
        "length": primer.length,
        "anneal_length": primer.anneal_length,
        "tm": primer.tm,
        "tm_full": primer.tm_full,
        "gc": primer.gc,
        "enzyme": primer.enzyme,
        "restriction_site": primer.restriction_site,
        "has_gc_clamp": primer.has_gc_clamp,
        "warnings": list(primer.warnings),
    }


def _end_payload(end) -> dict:
    return {"sequence": end.sequence, "strand": end.strand,
            "side": end.side, "kind": end.kind}


class PcrView(APIView):
    """POST /api/design/pcr/ — design a PCR and simulate its product.

    With no enzymes this is conventional PCR. Name a pair and each primer
    gains a tail carrying a site, the product is cut, and the insert that
    comes back is ready for `clone()`.
    """

    throttle_scope = "design"

    def post(self, request):
        serializer = PcrRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        try:
            result = design_pcr(
                data["template"],
                target_start=data["target_start"],
                target_end=data["target_end"],
                left_enzyme=data["left_enzyme"],
                right_enzyme=data["right_enzyme"],
                clamp=data["clamp"],
                keep_frame=data["keep_frame"],
                start_codon_mode=data["start_codon_mode"],
                name=data["name"],
                forward_primer=data["forward_primer"],
                reverse_primer=data["reverse_primer"],
            )
        except SequenceError as error:
            return _bad_request(error)

        payload = {
            "forward": _primer_payload(result.forward),
            "reverse": _primer_payload(result.reverse),
            "product": result.product,
            "product_length": result.product_length,
            "amplified_region": result.amplified_region,
            "template_start": result.template_start,
            "template_end": result.template_end,
            "annealing_temperature": result.annealing_temperature,
            "left_enzyme": result.left_enzyme,
            "right_enzyme": result.right_enzyme,
            "insert_orf_start": result.insert_orf_start,
            "problems": result.problems,
            "warnings": result.warnings,
            "is_clean": result.is_clean,
            "primer_source": "custom" if data["forward_primer"] is not None else "automatic",
            "digest": None,
            "gel": _gel_simulation(
                "Predicted PCR product",
                [
                    {
                        "name": "PCR",
                        "description": "Expected specific amplicon",
                        "bands": [{
                            "size_bp": result.product_length,
                            "label": f"{result.product_length:,} bp amplicon",
                        }],
                    },
                    {
                        "name": "NTC",
                        "description": "Expected no-template control",
                        "bands": [],
                    },
                ],
                [result.product_length],
            ),
            "preflight": pcr_preflight(
                result, keep_frame=data["keep_frame"],
            ).to_dict(),
            "provenance": _provenance("pcr", data, result.product),
        }
        if result.digest is not None:
            payload["digest"] = {
                "top": result.digest.top,
                "bottom": result.digest.bottom,
                "length": result.digest.length,
                "left_end": _end_payload(result.digest.left_end),
                "right_end": _end_payload(result.digest.right_end),
                "trimmed_left": result.digest.trimmed_left,
                "trimmed_right": result.digest.trimmed_right,
            }
        return Response(payload)


class FeatureDetectionView(APIView):
    """Read-only sequence annotation proposals for unsaved molecules."""

    def post(self, request):
        serializer = FeatureDetectionSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        try:
            matches = detect_common_features(
                data['sequence'], circular=data['circular'], existing=data['annotations'],
            )
        except SequenceError as error:
            return _bad_request(error)
        return Response({'matches': matches})
