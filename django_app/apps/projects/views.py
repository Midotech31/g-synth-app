"""Project CRUD — always scoped to the signed-in user."""
from django.http import HttpResponse
from django.utils import timezone
from rest_framework import viewsets
from rest_framework.decorators import action

from apps.projects.models import Project
from apps.projects.serializers import ProjectListSerializer, ProjectSerializer
from gsynth_engine.genbank import to_fasta, to_genbank
from gsynth_engine.provenance import build_provenance


class ProjectViewSet(viewsets.ModelViewSet):
    """
    Standard REST endpoints under /api/projects/:
        GET     /api/projects/       list current user's projects
        POST    /api/projects/       create
        GET     /api/projects/<id>/  retrieve
        PATCH   /api/projects/<id>/  partial update
        DELETE  /api/projects/<id>/  delete

    Every query is scoped to `request.user` — users can never see or touch
    another user's rows.
    """

    def get_queryset(self):
        # AnonymousUser is blocked earlier by DEFAULT_PERMISSION_CLASSES
        return Project.objects.filter(user=self.request.user)

    def get_serializer_class(self):
        return ProjectListSerializer if self.action == "list" else ProjectSerializer

    def perform_create(self, serializer):
        data = serializer.validated_data
        provenance = build_provenance(
            str(data.get("module", "general")),
            parameters={"name": data.get("name", ""), "sequence": data.get("sequence", "")},
            output_sequence=str(data.get("sequence", "")),
        )
        serializer.save(user=self.request.user, provenance=provenance)

    @action(detail=True, methods=["get"])
    def export(self, request, pk=None):
        """GET /api/projects/<id>/export/?filetype=genbank|fasta

        Saved work that cannot be taken out is not really saved. GenBank is
        the default because it is the only one of the two that carries the
        features, which is what the map is made of.
        """
        project = self.get_object()
        safe = (project.name or "sequence").replace(" ", "_")
        data = project.data or {}
        # Only the dedicated read-only field is authoritative. `data` is
        # user-editable JSON and must never be allowed to forge an audit hash.
        provenance = project.provenance or {}
        provenance_note = (
            f"G-Synth {provenance.get('engine_version', 'unknown')} · "
            f"output SHA-256 {provenance.get('output_sha256', 'unavailable')}"
        )

        # `format` is DRF's own renderer switch, so a value it does not
        # recognise 404s before this view ever runs.
        if request.query_params.get("filetype") == "fasta":
            description = " · ".join(filter(None, (project.notes, provenance_note)))
            body = to_fasta(project.sequence, name=safe, description=description)
            filename, content_type = f"{safe}.fasta", "text/plain; charset=utf-8"
        else:
            body = to_genbank(
                project.sequence,
                name=safe,
                description=project.notes or project.get_module_display(),
                features=data.get("annotations") or [],
                circular=str(data.get("topology", "")).lower() == "circular",
                date=timezone.now().strftime("%d-%b-%Y").upper(),
                comments=[provenance_note],
            )
            filename, content_type = f"{safe}.gb", "chemical/seq-na-genbank"

        response = HttpResponse(body, content_type=content_type)
        response["Content-Disposition"] = f'attachment; filename="{filename}"'
        return response
