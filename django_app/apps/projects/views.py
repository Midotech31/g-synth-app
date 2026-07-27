"""Project CRUD — always scoped to the signed-in user."""
from rest_framework import viewsets

from apps.projects.models import Project
from apps.projects.serializers import ProjectListSerializer, ProjectSerializer


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
        serializer.save(user=self.request.user)
