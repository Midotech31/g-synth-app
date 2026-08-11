"""Request and response shapes for the design endpoints.

The engine owns the biology; these serializers only validate what the user
sent and reshape what the engine returned. No design logic lives here — if
you find yourself computing a sequence in this file, it belongs in
`gsynth_engine`.
"""
from __future__ import annotations

from rest_framework import serializers

from gsynth_engine.constants import (
    ALL_ENZYMES,
    CLEAVAGE_SITES,
)
from gsynth_engine.merzoug import MAX_OVERHANG, MIN_OVERHANG
from gsynth_engine.vectors import CATALOGUE, DEFAULT_VECTOR

#: Every enzyme with verified cut geometry. The interface groups the
#: nineteen this lab keeps ahead of the rest; the API accepts all.
ENZYME_NAMES = sorted(ALL_ENZYMES)
CLEAVAGE_NAMES = sorted(CLEAVAGE_SITES)
VECTOR_KEYS = [spec.key for spec in CATALOGUE]


class DesignRequestSerializer(serializers.Serializer):
    """Common inputs for both SSD and assembly design."""

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
        """Translate the validated payload into engine arguments."""
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
    """SSD inputs plus the two knobs that control fragmentation."""

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
    """Opt in to persisting the design as a project."""

    save_as_project = serializers.BooleanField(default=False)


class SSDRequestSerializer(DesignRequestSerializer, SaveMixin):
    pass


class SaveableAssemblyRequestSerializer(AssemblyRequestSerializer, SaveMixin):
    pass


class VectorAnnotationSerializer(serializers.Serializer):
    """A feature carried over from a parsed vector file.

    Deliberately permissive about extra keys: the client sends back whatever
    the parser gave it, and dropping fields here would quietly lose colours
    and types on the way to the recombinant map.
    """

    name = serializers.CharField(max_length=200, required=False, default="")
    type = serializers.CharField(max_length=60, required=False, default="misc_feature")
    start = serializers.IntegerField(min_value=0)
    end = serializers.IntegerField(min_value=0)
    direction = serializers.IntegerField(required=False, default=0)
    color = serializers.CharField(max_length=32, required=False, default="")

    def validate(self, attrs):
        if attrs["end"] < attrs["start"]:
            raise serializers.ValidationError(
                {"end": "A feature cannot end before it starts."}
            )
        return attrs


class CloneRequestSerializer(AssemblyRequestSerializer, SaveMixin):
    """Design an insert and clone it into a vector, in one call.

    The insert is designed from the same inputs as /assembly/, so the two
    endpoints cannot disagree about what is being cloned.
    """

    #: Either name a catalogue vector, or supply a sequence — not neither.
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
    vector_is_circular = serializers.BooleanField(
        default=True,
        help_text="Only circular vectors can be cloned into — a linear one "
                  "cut twice leaves the backbone in two pieces.",
    )
    #: Design the full Merzoug assembly, or just the SSD duplex.
    fragment = serializers.BooleanField(
        default=True,
        help_text="False clones the SSD duplex directly, without fragmenting it.",
    )


def resolve_vector(data: dict) -> tuple[str, str, list[dict], object]:
    """Work out which vector to cut, from a key, a sequence, or both.

    Returns (sequence, name, annotations, spec). A supplied sequence always
    wins over the bundled one — a lab's own copy of a backbone is the one on
    their bench — but it is still checked against the catalogue entry, so a
    substitution is caught rather than silently cloned into.
    """
    from gsynth_engine import vectors as catalogue

    key = data.get("vector_key") or ""
    supplied = (data.get("vector") or "").strip()
    spec = catalogue.get(key) if key else None

    if supplied:
        sequence = supplied
        # An unnamed sequence may still be recognisable.
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
    return sequence, name, annotations, spec


class OptimiseRequestSerializer(serializers.Serializer):
    """Rewrite a gene for a host, keeping the protein identical."""

    sequence = serializers.CharField(
        max_length=200_000,
        help_text="A coding sequence, or a protein when is_protein is set.",
    )
    is_protein = serializers.BooleanField(default=False)
    keep_stop = serializers.BooleanField(
        default=True,
        help_text="Turn off for an insert destined for a C-terminal vector "
                  "tag, where a stop codon would silently remove the tag.",
    )
    #: The cloning pair, so the optimiser removes sites that would break it.
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
    #: Reference genes to measure the usage table against, when the CAI matters.
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
    """How much of each thing goes in the tube."""

    vector_length = serializers.IntegerField(min_value=1, max_value=1_000_000)
    insert_length = serializers.IntegerField(min_value=1, max_value=1_000_000)
    vector_ng = serializers.FloatField(min_value=0.01, max_value=100_000, default=50.0)
    ends = serializers.ChoiceField(choices=["5'", "3'", "blunt"], default="5'")
    total_volume_uL = serializers.FloatField(min_value=1, max_value=1000, default=20.0)
    #: Left empty, the recommendation for the kind of ends is used.
    ratios = serializers.ListField(
        child=serializers.FloatField(min_value=0.01, max_value=100),
        required=False, default=list, max_length=8,
    )


class PrimerRequestSerializer(serializers.Serializer):
    """Sequencing primers that read across a region."""

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
    """Compare sequencing reads to the design."""

    design = serializers.CharField(max_length=2_000_000)
    #: Named so a report can say which trace disagreed.
    reads = serializers.DictField(child=serializers.CharField(max_length=100_000))
    circular = serializers.BooleanField(default=True)
    trim = serializers.IntegerField(min_value=0, max_value=200, default=30)
    coding_start = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    coding_end = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    region_start = serializers.IntegerField(min_value=0, required=False, allow_null=True)
    region_end = serializers.IntegerField(min_value=0, required=False, allow_null=True)

    #: DictField has no max_length, so the cap is enforced here.
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
    """Sanger trace files, compared against a design.

    A trace carries what the letters cannot: how much to believe each base.
    Files arrive as multipart because an .ab1 is binary — base64 in JSON
    would inflate a 400 kB trace to 550 kB for no gain.
    """

    design = serializers.CharField(max_length=2_000_000)
    #: Bounded, because an .ab1 is a fixed shape and anything far larger is
    #: either the wrong file or someone probing for a memory limit.
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

    #: A 1.2 kb read is about 400 kB. Ten times that is not a trace.
    MAX_TRACE_BYTES = 4 * 1024 * 1024

    def validate_traces(self, files):
        for upload in files:
            if upload.size > self.MAX_TRACE_BYTES:
                raise serializers.ValidationError(
                    f"{upload.name} is {upload.size / 1e6:.1f} MB. A Sanger "
                    f"trace is normally under 0.5 MB — check it is an .ab1 "
                    f"and not an archive."
                )
        return files


class AlignRequestSerializer(serializers.Serializer):
    """Compare two sequences that are not assumed to be the same thing."""

    first = serializers.CharField(max_length=200_000)
    second = serializers.CharField(max_length=200_000)
    mode = serializers.ChoiceField(
        choices=["global", "local", "semi-global"], default="global",
    )
    is_protein = serializers.BooleanField(default=False)
    #: A gene cloned the other way round is not a different gene.
    try_reverse = serializers.BooleanField(default=True)
    match = serializers.IntegerField(default=5, min_value=1, max_value=20)
    mismatch = serializers.IntegerField(default=-4, min_value=-20, max_value=0)
    gap_open = serializers.IntegerField(default=10, min_value=0, max_value=60)
    gap_extend = serializers.IntegerField(default=1, min_value=0, max_value=20)
