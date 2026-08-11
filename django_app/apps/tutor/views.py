"""The study assistant — a chat proxy in front of a local Ollama model.

Distinct from every other view in this codebase: it does not call
`gsynth_engine`, because there is no deterministic biology to compute here.
What it returns is a language model's answer, not a verified result — the
request is bounded the same way every other endpoint's is, but nothing
about the response is checked the way a design or a clone is.
"""
from __future__ import annotations

from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.tutor.ollama import OllamaError, ask
from apps.tutor.serializers import TutorRequestSerializer


class TutorView(APIView):
    """POST /api/tutor/ask/ — ask the study assistant a question.

    503, not 400 or 500, when Ollama cannot be reached: the request itself
    was fine, the model just is not there to answer it — the same distinction
    a database being down gets over a query being malformed.
    """

    throttle_scope = "tutor"

    def post(self, request):
        serializer = TutorRequestSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        history = [{"role": turn["role"], "content": turn["content"]} for turn in data["history"]]
        try:
            answer = ask(data["question"], history)
        except OllamaError as error:
            return Response({"detail": str(error)}, status=status.HTTP_503_SERVICE_UNAVAILABLE)

        return Response({"answer": answer})
