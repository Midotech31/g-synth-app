"""Client for an optional local Ollama service."""
from __future__ import annotations

import json
import urllib.error
import urllib.request

from django.conf import settings


class OllamaError(Exception):
    """Ollama could not be reached or refused the request."""


def ask(question: str, history: list[dict]) -> str:
    """Send one bounded conversation turn to the configured model."""
    messages = [*history, {"role": "user", "content": question}]

    body = json.dumps({
        "model": settings.OLLAMA_MODEL,
        "messages": messages,
        "stream": False,
    }).encode("utf-8")

    request = urllib.request.Request(
        f"{settings.OLLAMA_BASE_URL}/api/chat",
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=settings.OLLAMA_TIMEOUT_SECONDS) as response:
            payload = json.loads(response.read())
    except urllib.error.HTTPError as error:
        try:
            detail = json.loads(error.read()).get("error", "")
        except (ValueError, AttributeError):
            detail = ""
        hint = (f" Try `ollama pull {settings.OLLAMA_MODEL}` first."
                if "not found" in detail else "")
        raise OllamaError(
            f"Ollama refused the request{f': {detail}' if detail else '.'}{hint}"
        ) from error
    except (urllib.error.URLError, TimeoutError) as error:
        raise OllamaError(
            f"Could not reach Ollama at {settings.OLLAMA_BASE_URL} within "
            f"{settings.OLLAMA_TIMEOUT_SECONDS}s. Make sure `ollama serve` is "
            f"running and that {settings.OLLAMA_MODEL!r} has been pulled "
            f"(`ollama pull {settings.OLLAMA_MODEL}`)."
        ) from error
    except json.JSONDecodeError as error:
        raise OllamaError(
            "Ollama returned something that was not JSON. Check that "
            "OLLAMA_BASE_URL points at Ollama's own API, not something "
            "else on that port."
        ) from error

    message = payload.get("message", {}).get("content", "").strip()
    if not message:
        raise OllamaError("Ollama returned an empty answer. Try asking again.")
    return message
