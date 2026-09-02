"""Tutor routes — mounted under /api/tutor/ in config.urls."""
from django.urls import path

from apps.tutor.views import TutorStatusView, TutorView

urlpatterns = [
    path("status/", TutorStatusView.as_view(), name="tutor-status"),
    path("ask/", TutorView.as_view(), name="tutor-ask"),
]
