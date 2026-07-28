"""Request and response shapes for the design endpoints.

The engine owns the biology; these serializers only validate what the user
sent and reshape what the engine returned. No design logic lives here — if
you find yourself computing a sequence in this file, it belongs in
`gsynth_engine`.
"""
from __future__ import annotations

from gsynth_engine.constants import (
    CLEAVAGE_SITES,
    RESTRICTION_ENZYMES,
)
from gsynth_engine.merzoug import MAX_OVERHANG, MIN_OVERHANG
from gsynth_engine.vectors import CATALOGUE, DEFAULT_VECTOR
from rest_framework import serializers

ENZYME_NAMES = sorted(RESTRICTION_ENZYMES)
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
