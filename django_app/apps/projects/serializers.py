from rest_framework import serializers

from apps.projects.models import Project


class ProjectSerializer(serializers.ModelSerializer):
    """Full project representation for retrieve / update."""

    class Meta:
        model = Project
        fields = ("id", "name", "module", "sequence", "notes", "data",
                  "created_at", "updated_at")
        read_only_fields = ("id", "created_at", "updated_at")


class ProjectListSerializer(serializers.ModelSerializer):
    """Compact projection for the list endpoint — omits `data` and `sequence`.

    Keeps the list payload small even when a user has hundreds of projects
    with heavy JSON payloads. Callers fetch full details via retrieve.
    """

    class Meta:
        model = Project
        fields = ("id", "name", "module", "updated_at")
        read_only_fields = fields
