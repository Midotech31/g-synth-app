"""The study assistant — a chat proxy in front of a local Ollama model.

Distinct from every other view in this codebase: it does not call
`gsynth_engine`, because there is no deterministic biology to compute here.
What it returns is a language model's answer, not a verified result — the
request is bounded the same way every other endpoint's is, but nothing
about the response is checked the way a design or a clone is.
"""
from __future__ import annotations

from django.conf import settings
from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.tutor.ollama import OllamaError, ask
from apps.tutor.serializers import TutorRequestSerializer

TUTOR_NOTICE = (
    "This optional AI study assistant is not a validated scientific result. "
    "Do not enter confidential or unpublished sequences unless you administer "
    "and trust the configured Ollama service."
)


class TutorStatusView(APIView):
    """GET /api/tutor/status/ — disclose availability before data entry."""

    def get(self, request):
        enabled = bool(settings.TUTOR_ENABLED)
        return Response({
            "enabled": enabled,
            "notice": TUTOR_NOTICE,
            "disabled_reason": (
                "Learn is disabled on this deployment. The deterministic "
                "G-Synth design and verification tools remain available."
                if not enabled else ""
            ),
        })


class TutorView(APIView):
    """POST /api/tutor/ask/ — ask the study assistant a question.

    503, not 400 or 500, when Ollama cannot be reached: the request itself
    was fine, the model just is not there to answer it — the same distinction
    a database being down gets over a query being malformed.
    """

    throttle_scope = "tutor"

    def post(self, request):
        if not settings.TUTOR_ENABLED:
            return Response(
                {
                    "code": "tutor_disabled",
                    "detail": (
                        "Learn is disabled on this deployment. The deterministic "
                        "G-Synth design and verification tools remain available."
                    ),
                },
                status=status.HTTP_503_SERVICE_UNAVAILABLE,
            )
        serializer = TutorRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        history = [{"role": turn["role"], "content": turn["content"]} for turn in data["history"]]
        try:
            answer = ask(data["question"], history)
        except OllamaError as error:
            return Response({"detail": str(error)}, status=status.HTTP_503_SERVICE_UNAVAILABLE)

        return Response({"answer": answer})
