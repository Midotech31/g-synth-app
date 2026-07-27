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
from rest_framework import serializers

ENZYME_NAMES = sorted(RESTRICTION_ENZYMES)
CLEAVAGE_NAMES = sorted(CLEAVAGE_SITES)


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
