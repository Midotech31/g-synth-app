from __future__ import annotations

from rest_framework import serializers

from gsynth_engine.codon import DEFAULT_HOST, TABLES
from gsynth_engine.constants import (
    ALL_ENZYMES,
    CLEAVAGE_SITES,
)
from gsynth_engine.esd import MAX_OVERHANG, MIN_OVERHANG
from gsynth_engine.pcr import DEFAULT_CLAMP
from gsynth_engine.vectors import CATALOGUE, DEFAULT_VECTOR

ENZYME_NAMES = sorted(ALL_ENZYMES)
CLEAVAGE_NAMES = sorted(CLEAVAGE_SITES)
VECTOR_KEYS = [spec.key for spec in CATALOGUE]


class DesignRequestSerializer(serializers.Serializer):


    sequence = serializers.CharField(
        max_length=200_000,
        help_text="The insert, A/C/G/T. FASTA headers and whitespace are stripped.",
    )
    name = serializers.CharField(
        max_length=120, required=False, default="construct",
        help_text="Used to name the oligos in the order sheet.",
    )
    left_enzyme = serializers.ChoiceField(choices=ENZYME_NAMES, default="NdeI")
    right_enzyme = serializers.ChoiceField(choices=ENZYME_NAMES, default="XhoI")
    is_coding = serializers.BooleanField(
        default=False,
        help_text="True when the insert already carries its own ATG.",
    )
    remove_stop = serializers.BooleanField(
        default=False,
        help_text="Coding inserts only: truncate at the first in-frame stop.",
    )
    cleavage_site = serializers.ChoiceField(
        choices=CLEAVAGE_NAMES, required=False, allow_null=True,
        allow_blank=True, default="Thrombin",
    )
    include_his_tag = serializers.BooleanField(default=True)
    include_linkers = serializers.BooleanField(default=True)

    def validate(self, attrs):
        if attrs["left_enzyme"] == attrs["right_enzyme"]:
            raise serializers.ValidationError({
                "right_enzyme": "The two enzymes must differ, otherwise the "
                                "insert could ligate in either orientation.",
            })
        return attrs

    @property
    def engine_kwargs(self) -> dict:

        data = self.validated_data
        return {
            "enzyme_pair": f"{data['left_enzyme']} / {data['right_enzyme']}",
            "is_coding": data["is_coding"],
            "remove_stop": data["remove_stop"],
            "cleavage_site": data.get("cleavage_site") or None,
            "include_his_tag": data["include_his_tag"],
            "include_linkers": data["include_linkers"],
        }


class AssemblyRequestSerializer(DesignRequestSerializer):


    target_oligo_length = serializers.IntegerField(
        min_value=20, max_value=300, default=90,
        help_text="How long each synthesised oligo should be.",
    )
    overhang_length = serializers.IntegerField(
        min_value=MIN_OVERHANG, max_value=MAX_OVERHANG, default=4,
        help_text=f"Junction overhang, {MIN_OVERHANG}–{MAX_OVERHANG} nt.",
    )

    @property
    def engine_kwargs(self) -> dict:
        kwargs = super().engine_kwargs
        kwargs["target_oligo_length"] = self.validated_data["target_oligo_length"]
        kwargs["overhang_length"] = self.validated_data["overhang_length"]
        return kwargs


class SaveMixin(serializers.Serializer):


    save_as_project = serializers.BooleanField(default=False)


class SSDRequestSerializer(DesignRequestSerializer, SaveMixin):
    pass


class SaveableAssemblyRequestSerializer(AssemblyRequestSerializer, SaveMixin):
    pass


class VectorAnnotationSerializer(serializers.Serializer):


    name = serializers.CharField(max_length=200, required=False, default="")
    type = serializers.CharField(max_length=60, required=False, default="misc_feature")
    start = serializers.IntegerField(min_value=0)
    end = serializers.IntegerField(min_value=0)
    direction = serializers.IntegerField(required=False, default=0)
    color = serializers.CharField(max_length=32, required=False, default="")
    translation_start = serializers.IntegerField(min_value=0, required=False)
    translation_end = serializers.IntegerField(min_value=0, required=False)
    truncated = serializers.BooleanField(required=False)
    inferred = serializers.BooleanField(required=False)
    basis = serializers.CharField(max_length=4000, required=False, allow_blank=True)
    regulatory_class = serializers.CharField(max_length=80, required=False)

    def validate(self, attrs):
        if attrs["end"] < attrs["start"]:
            raise serializers.ValidationError(
                {"end": "A feature cannot end before it starts."}
            )
        translation_start = attrs.get("translation_start")
        translation_end = attrs.get("translation_end")
        if (translation_start is None) != (translation_end is None):
            raise serializers.ValidationError(
                "Translation start and end must be supplied together."
            )
        if translation_start is not None and (
            translation_start < attrs["start"]
            or translation_end > attrs["end"]
            or translation_end <= translation_start
        ):
            raise serializers.ValidationError(
                "Translation coordinates must define a non-empty span inside the feature."
            )
        return attrs


class CloneRequestSerializer(AssemblyRequestSerializer, SaveMixin):


    vector_key = serializers.ChoiceField(
        choices=VECTOR_KEYS, required=False, allow_blank=True,
        default=DEFAULT_VECTOR.key,
        help_text="A catalogue vector. Its bundled sequence is used unless one "
                  "is supplied, and either way the supplied sequence is checked "
                  "against the entry.",
    )
    vector = serializers.CharField(
        max_length=2_000_000, required=False, allow_blank=True, default="",
        help_text="The vector's full sequence. Circular; the origin can be "
                  "anywhere. Optional when vector_key names a bundled vector.",
    )
    vector_name = serializers.CharField(max_length=200, required=False, default="")
    vector_annotations = VectorAnnotationSerializer(many=True, required=False)
    product_annotations = VectorAnnotationSerializer(many=True, required=False)
    insert_annotations = VectorAnnotationSerializer(many=True, required=False, max_length=5000)
    vector_is_circular = serializers.BooleanField(
        default=True,
        help_text="Only circular vectors can be cloned into — a linear one "
                  "cut twice leaves the backbone in two pieces.",
    )

    fragment = serializers.BooleanField(
        default=True,
        help_text="False clones the SSD duplex directly, without fragmenting it.",
    )

    pre_digested = serializers.BooleanField(default=False)
    orf_start = serializers.IntegerField(
        min_value=0,
        required=False,
        allow_null=True,
        default=None,
        help_text="Known zero-based translation start in a supplied digested insert.",
    )
    insert_reverse = serializers.CharField(
        max_length=200_000, required=False, allow_blank=True, default="",
    )

    def validate(self, attrs):
        attrs = super().validate(attrs)
        if attrs.get("insert_annotations"):
            if not attrs.get("pre_digested"):
                raise serializers.ValidationError({"insert_annotations": "Insert annotations require a supplied duplex."})
            for annotation in attrs["insert_annotations"]:
                if not 0 <= annotation["start"] < annotation["end"] <= len(attrs["sequence"]):
                    raise serializers.ValidationError({"insert_annotations": "Feature coordinates must lie within the supplied top strand."})
        if attrs.get("pre_digested") and not attrs.get("insert_reverse"):
            raise serializers.ValidationError({
                "insert_reverse": "Supply the cut fragment's bottom strand as "
                                  "well. With only one strand the overhangs "
                                  "cannot be measured, and an insert cut for a "
                                  "different pair would look correct."
            })
        if (
            attrs.get("pre_digested")
            and attrs.get("orf_start") is not None
            and attrs["orf_start"] + 3 > len(attrs["sequence"])
        ):
            raise serializers.ValidationError({
                "orf_start": "The translation start must identify a complete codon inside the insert."
            })
        for annotation in attrs.get("product_annotations") or []:
            if annotation["end"] == annotation["start"]:
                raise serializers.ValidationError({
                    "product_annotations": "A product feature must cover at least one base."
                })
        return attrs


def resolve_vector(data: dict) -> tuple[str, str, list[dict], object]:

    from gsynth_engine import vectors as catalogue
    from gsynth_engine.annotations import detect_common_features

    key = data.get("vector_key") or ""
    supplied = (data.get("vector") or "").strip()
    spec = catalogue.get(key) if key else None

    if supplied:
        sequence = supplied

        spec = spec or catalogue.identify(supplied)
    else:
        record = catalogue.sequence_of(spec.key) if spec else None
        if record is None:
            raise serializers.ValidationError({
                "vector": "No vector sequence. Choose a vector whose sequence "
                          "ships with G-Synth, or paste or import your own.",
            })
        sequence = record["sequence"]
        if not data.get("vector_annotations"):
            data["vector_annotations"] = record["annotations"]

    name = data.get("vector_name") or (spec.name if spec else "vector")
    annotations = [dict(f) for f in (data.get("vector_annotations") or [])]
    detected = detect_common_features(sequence, circular=True, existing=annotations)
    for match in detected:
        annotation = dict(match["annotation"])
        if annotation["type"] not in {"promoter", "RBS", "terminator", "protein_bind"}:
            continue
        annotation["inferred"] = True
        annotation["basis"] = match["basis"]
        annotations.append(annotation)
    annotations.sort(key=lambda feature: (feature["start"], feature["end"], feature["name"]))
    return sequence, name, annotations, spec


class FeatureDetectionSerializer(serializers.Serializer):
    sequence = serializers.CharField(max_length=200_000)
    circular = serializers.BooleanField(default=False)
    annotations = VectorAnnotationSerializer(many=True, required=False, default=list, max_length=5000)


class OptimiseRequestSerializer(serializers.Serializer):


    sequence = serializers.CharField(
        max_length=200_000,
        help_text="A coding sequence, or a protein when is_protein is set.",
    )
    host = serializers.ChoiceField(choices=tuple(TABLES), default=DEFAULT_HOST)
    is_protein = serializers.BooleanField(default=False)
    protein_context = serializers.ChoiceField(
        choices=("auto", "mature_peptide", "complete_orf"),
        default="auto",
        help_text="How an N-terminal methionine is handled during peptide-to-DNA conversion.",
    )
    keep_stop = serializers.BooleanField(
        default=True,
        help_text="Turn off for an insert destined for a C-terminal vector "
                  "tag, where a stop codon would silently remove the tag.",
    )

    avoid_enzymes = serializers.ListField(
        child=serializers.ChoiceField(choices=ENZYME_NAMES),
        required=False, default=list, max_length=12,
    )
    avoid_motifs = serializers.ListField(
        child=serializers.RegexField(r"^[ACGTacgt]{2,40}$"),
        required=False, default=list, max_length=20,
    )
    max_homopolymer = serializers.IntegerField(min_value=3, max_value=12, default=5)
    gc_min = serializers.FloatField(min_value=10.0, max_value=90.0, default=30.0)
    gc_max = serializers.FloatField(min_value=10.0, max_value=90.0, default=70.0)
    gc_window = serializers.IntegerField(min_value=20, max_value=200, default=50)
    max_repeat = serializers.IntegerField(min_value=8, max_value=40, default=15)
    avoid_rare = serializers.BooleanField(default=True)

    reference_genes = serializers.ListField(
        child=serializers.CharField(max_length=100_000),
        required=False, default=list, max_length=200,
    )

    def validate(self, attrs):
        if attrs["gc_min"] >= attrs["gc_max"]:
            raise serializers.ValidationError({
                "gc_max": "The upper GC bound must be above the lower one.",
            })
        return attrs


class LigationRequestSerializer(serializers.Serializer):


    vector_length = serializers.IntegerField(min_value=1, max_value=1_000_000)
    insert_length = serializers.IntegerField(min_value=1, max_value=1_000_000)
    vector_ng = serializers.FloatField(min_value=0.01, max_value=100_000, default=50.0)
    ends = serializers.ChoiceField(choices=["5'", "3'", "blunt"], default="5'")
    total_volume_uL = serializers.FloatField(min_value=1, max_value=1000, default=20.0)

    ratios = serializers.ListField(
        child=serializers.FloatField(min_value=0.01, max_value=100),
        required=False, default=list, max_length=8,
    )


class PrimerRequestSerializer(serializers.Serializer):


    template = serializers.CharField(max_length=2_000_000)
    target_start = serializers.IntegerField(min_value=0)
    target_end = serializers.IntegerField(min_value=1)
    circular = serializers.BooleanField(default=True)
    name = serializers.CharField(max_length=60, required=False, default="seq")
    tm_min = serializers.FloatField(min_value=30, max_value=85, default=50.0)
    tm_max = serializers.FloatField(min_value=30, max_value=90, default=65.0)
    margin = serializers.IntegerField(min_value=50, max_value=2000, default=80)
    read_length = serializers.IntegerField(min_value=100, max_value=1500, default=700)

    def validate(self, attrs):
        if attrs["target_end"] <= attrs["target_start"]:
            raise serializers.ValidationError(
                {"target_end": "The region must end after it starts."}
            )
        if attrs["tm_max"] <= attrs["tm_min"]:
            raise serializers.ValidationError(
                {"tm_max": "The upper Tm bound must be above the lower one."}
            )
        return attrs


class VerifyRequestSerializer(serializers.Serializer):


    design = serializers.CharField(max_length=2_000_000)

    reads = serializers.DictField(child=serializers.CharField(max_length=100_000))
    circular = serializers.BooleanField(default=True)
    trim = serializers.IntegerField(min_value=0, max_value=200, default=30)
    coding_start = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    coding_end = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    region_start = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    region_end = serializers.IntegerField(min_value=0, required=False, allow_null=True)


    MAX_READS = 64

    def validate(self, attrs):
        if not attrs["reads"]:
            raise serializers.ValidationError({"reads": "Send at least one read."})
        if len(attrs["reads"]) > self.MAX_READS:
            raise serializers.ValidationError(
                {"reads": f"At most {self.MAX_READS} reads at a time."}
            )
        return attrs


class TraceUploadSerializer(serializers.Serializer):


    design = serializers.CharField(max_length=2_000_000)


    traces = serializers.ListField(
        child=serializers.FileField(max_length=255, allow_empty_file=False),
        min_length=1, max_length=32,
    )
    circular = serializers.BooleanField(default=True)
    trim_quality = serializers.IntegerField(min_value=0, max_value=60, default=13)
    coding_start = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    coding_end = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    region_start = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    region_end = serializers.IntegerField(min_value=0, required=False, allow_null=True)


    MAX_TRACE_BYTES = 4 * 1024 * 1024

    def validate_traces(self, files):
        for upload in files:
            if upload.size > self.MAX_TRACE_BYTES:
                raise serializers.ValidationError(
                    f"{upload.name} is {upload.size / 1e6:.1f} MB. A Sanger "
                    f"trace is normally under 0.5 MB — check it is an .ab1 or .scf "
                    f"and not an archive."
                )
        return files


class AlignRequestSerializer(serializers.Serializer):


    first = serializers.CharField(max_length=200_000)
    second = serializers.CharField(max_length=200_000)
    mode = serializers.ChoiceField(
        choices=["global", "local", "semi-global"], default="global",
    )
    is_protein = serializers.BooleanField(default=False)

    try_reverse = serializers.BooleanField(default=True)
    match = serializers.IntegerField(default=5, min_value=1, max_value=20)
    mismatch = serializers.IntegerField(default=-4, min_value=-20, max_value=0)
    gap_open = serializers.IntegerField(default=10, min_value=0, max_value=60)
    gap_extend = serializers.IntegerField(default=1, min_value=0, max_value=20)


class HybridizationRequestSerializer(serializers.Serializer):


    first = serializers.CharField(max_length=200_000)
    second = serializers.CharField(max_length=200_000)
    analysis_temperature_c = serializers.FloatField(
        default=25.0, min_value=-100.0, max_value=150.0,
    )
    oligo_nM = serializers.FloatField(
        default=50_000.0, min_value=0.001, max_value=100_000_000.0,
        help_text="Total strand concentration used by the Tm model.",
    )
    na_mM = serializers.FloatField(default=50.0, min_value=0.0, max_value=5_000.0)
    mg_mM = serializers.FloatField(default=0.0, min_value=0.0, max_value=1_000.0)
    dntp_mM = serializers.FloatField(default=0.0, min_value=0.0, max_value=100.0)


class PcrRequestSerializer(serializers.Serializer):


    template = serializers.CharField(max_length=200_000)
    target_start = serializers.IntegerField(min_value=0, default=0)
    target_end = serializers.IntegerField(min_value=1, required=False, allow_null=True,
                                          default=None)
    left_enzyme = serializers.ChoiceField(choices=ENZYME_NAMES, required=False,
                                          allow_null=True, default=None)
    right_enzyme = serializers.ChoiceField(choices=ENZYME_NAMES, required=False,
                                           allow_null=True, default=None)

    clamp = serializers.IntegerField(min_value=0, max_value=20, default=DEFAULT_CLAMP)
    keep_frame = serializers.BooleanField(default=False)
    start_codon_mode = serializers.ChoiceField(
        choices=("use_site", "keep_both"), default="use_site",
        help_text=(
            "When the left site supplies ATG, use that start codon by default "
            "or deliberately retain the template ATG as well."
        ),
    )
    name = serializers.CharField(max_length=60, required=False, default="product")
    forward_primer = serializers.CharField(
        max_length=300, required=False, allow_null=True, default=None,
    )
    reverse_primer = serializers.CharField(
        max_length=300, required=False, allow_null=True, default=None,
    )

    def validate(self, attrs):
        end = attrs.get("target_end")
        if end is not None and end <= attrs["target_start"]:
            raise serializers.ValidationError(
                {"target_end": "The region must end after it starts."}
            )
        if (attrs.get("left_enzyme") is None) != (attrs.get("right_enzyme") is None):
            raise serializers.ValidationError(
                {"right_enzyme": "Choose an enzyme for both ends, or for neither."}
            )
        if (
            attrs.get("left_enzyme") is not None
            and attrs.get("left_enzyme") == attrs.get("right_enzyme")
        ):
            raise serializers.ValidationError({
                "right_enzyme": "The two enzymes must differ, otherwise the "
                                "insert could ligate in either orientation."
            })
        if (attrs.get("forward_primer") is None) != (attrs.get("reverse_primer") is None):
            raise serializers.ValidationError({
                "reverse_primer": "Enter both custom primer sequences, or neither."
            })
        return attrs
