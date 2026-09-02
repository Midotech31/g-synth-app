"""Optional private study-assistant endpoints."""
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
    """Report study-assistant availability and scientific limits."""

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
    """Submit a bounded question to the configured private service."""

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
