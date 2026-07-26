"""Top-level URL routes."""
from django.contrib import admin
from django.http import JsonResponse
from django.urls import include, path


def health(_request):
    """Liveness probe used by Docker HEALTHCHECK, Render, load balancers."""
    return JsonResponse({"status": "ok", "service": "gsynth-api"})


urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/health/", health, name="health"),
    path("api/auth/", include("apps.accounts.urls")),
    path("api/projects/", include("apps.projects.urls")),
]
