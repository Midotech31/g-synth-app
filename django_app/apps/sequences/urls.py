from django.urls import path

from apps.sequences.views import ParseSequenceFileView

urlpatterns = [
    path("parse/", ParseSequenceFileView.as_view(), name="sequence-parse"),
]
