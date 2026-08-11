"""A thin client for a local Ollama server.

Deliberately not a third-party HTTP client: this is one JSON POST with a
timeout, and `urllib` in the standard library does that without adding a
dependency whose only job in this codebase would be this file.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.request

from django.conf import settings

#: Frames the model as this app's study assistant rather than a bare chat
#: endpoint, and says what it should decline — a question with no answer is
#: better than a confident one about something else entirely.
SYSTEM_PROMPT = (
    "You are the study assistant built into G-Synth, a gene synthesis and "
    "cloning tool. Answer questions about molecular biology, gene synthesis, "
    "cloning, and genetic engineering clearly and accurately, at the level "
    "of a graduate student or research technician. Where it genuinely helps, "
    "relate the answer to G-Synth's own workflow: Optimise rewrites codon "
    "usage for a host without changing the protein; Design builds a cassette "
    "and cuts it into oligo pairs for Merzoug assembly, hybridised with "
    "sticky-end overhangs and ligated in order with no PCR at any step; "
    "Clone cuts a vector and ligates the insert in; Check works out ligation "
    "ratios, sequencing primers, and reads .ab1 traces back against the "
    "design; Compare aligns two sequences that are not assumed to be the "
    "same thing. Do not force this connection where a question does not "
    "call for it. If a question falls outside molecular biology and genetic "
    "engineering, say so briefly rather than answering it anyway. Keep "
    "answers concise — a few short paragraphs unless the question asks for "
    "more."
)


class OllamaError(Exception):
    """Ollama could not be reached or refused the request.

    The message is written for the person reading it on screen, the same
    convention as `SequenceError` — what is wrong, and what to do about it.
    """


def ask(question: str, history: list[dict]) -> str:
    """One turn of conversation with the configured Ollama model.

    `history` is prior turns, oldest first, each `{"role", "content"}` —
    bounding its size is the serialiser's job, not this function's.
    """
    messages = [{"role": "system", "content": SYSTEM_PROMPT}, *history,
                {"role": "user", "content": question}]

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
