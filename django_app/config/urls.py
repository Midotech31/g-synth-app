from django.contrib import admin
from django.http import JsonResponse
from django.urls import include, path


def health(_request):

    return JsonResponse({"status": "ok", "service": "gsynth-api"})


urlpatterns = [
    path("admin/", admin.site.urls),
    path("api/health/", health, name="health"),
    path("api/auth/", include("apps.accounts.urls")),
    path("api/projects/", include("apps.projects.urls")),
    path("api/sequences/", include("apps.sequences.urls")),
    path("api/design/", include("apps.design.urls")),
    path("api/tutor/", include("apps.tutor.urls")),
]
