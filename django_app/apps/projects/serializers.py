from types import SimpleNamespace

from rest_framework import serializers

from apps.projects.models import Project


class AnnotationSerializer(serializers.Serializer):
    """One user-editable sequence feature in engine coordinates."""

    name = serializers.CharField(max_length=200, trim_whitespace=True)
    type = serializers.CharField(max_length=64, default="misc_feature")
    start = serializers.IntegerField(min_value=0)
    end = serializers.IntegerField(min_value=1)
    direction = serializers.ChoiceField(choices=(-1, 0, 1), default=1)
    color = serializers.RegexField(
        r"^#[0-9A-Fa-f]{6}$", default="#0E6E77",
        error_messages={"invalid": "Use a six-digit hexadecimal colour such as #0E6E77."},
    )
    truncated = serializers.BooleanField(required=False)
    inferred = serializers.BooleanField(required=False)
    basis = serializers.CharField(max_length=4000, required=False, allow_blank=True)
    regulatory_class = serializers.CharField(max_length=80, required=False)
    translation_start = serializers.IntegerField(min_value=0, required=False)
    translation_end = serializers.IntegerField(min_value=1, required=False)

    def validate(self, attrs):
        project = self.context["project"]
        sequence_length = len(project.sequence or "")
        start, end = attrs["start"], attrs["end"]
        circular = str((project.data or {}).get("topology", "")).lower() == "circular"

        if not sequence_length:
            raise serializers.ValidationError("Annotations require a non-empty project sequence.")
        if start >= sequence_length:
            raise serializers.ValidationError({"start": "Start must fall inside the sequence."})
        if end <= start:
            raise serializers.ValidationError({"end": "End must be after start."})
        if end - start > sequence_length:
            raise serializers.ValidationError({"end": "A feature cannot be longer than the sequence."})
        if not circular and end > sequence_length:
            raise serializers.ValidationError({"end": "A linear feature cannot cross the sequence end."})

        translation_start = attrs.get("translation_start")
        translation_end = attrs.get("translation_end")
        if (translation_start is None) != (translation_end is None):
            raise serializers.ValidationError(
                "Translation start and end must either both be supplied or both be omitted."
            )
        if translation_start is not None and not (
            start <= translation_start < translation_end <= end
        ):
            raise serializers.ValidationError(
                "Translation coordinates must lie inside the feature."
            )
        return attrs


class ProjectAnnotationsSerializer(serializers.Serializer):
    """Validated replacement for only the annotations inside project data."""

    annotations = AnnotationSerializer(many=True)
    expected_updated_at = serializers.DateTimeField(required=False)

    def validate_annotations(self, value):
        if len(value) > 5000:
            raise serializers.ValidationError("A project can contain at most 5,000 features.")
        return value


class ProjectSerializer(serializers.ModelSerializer):
    """Full project representation for retrieve / update."""

    def validate(self, attrs):
        if self.instance and 'sequence' in attrs and attrs['sequence'] != self.instance.sequence:
            raise serializers.ValidationError({'sequence': (
                'Save the changed sequence as a new project. The original sequence '
                'must remain bound to its saved results and provenance.'
            )})
        data = attrs.get('data', self.instance.data if self.instance else {})
        if not isinstance(data, dict):
            raise serializers.ValidationError({'data': 'Project data must be an object.'})
        if 'data' in attrs and 'annotations' in data:
            holder = SimpleNamespace(
                sequence=attrs.get('sequence', self.instance.sequence if self.instance else ''),
                data=data,
            )
            features = ProjectAnnotationsSerializer(
                data={'annotations': data['annotations']}, context={'project': holder},
            )
            features.is_valid(raise_exception=True)
            attrs['data'] = {**data, 'annotations': features.validated_data['annotations']}
        return attrs

    class Meta:
        model = Project
        fields = ("id", "name", "module", "sequence", "notes", "data", "provenance",
                  "created_at", "updated_at")
        read_only_fields = ("id", "provenance", "created_at", "updated_at")


class ProjectListSerializer(serializers.ModelSerializer):
    """Compact projection for the list endpoint — omits `data` and `sequence`.

    Keeps the list payload small even when a user has hundreds of projects
    with heavy JSON payloads. Callers fetch full details via retrieve.
    """

    class Meta:
        model = Project
        fields = ("id", "name", "module", "updated_at")
        read_only_fields = fields
