"""Tutor routes — mounted under /api/tutor/ in config.urls."""
from django.urls import path

from apps.tutor.views import TutorView

urlpatterns = [
    path("ask/", TutorView.as_view(), name="tutor-ask"),
]
