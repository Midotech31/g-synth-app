"""Tests for the study assistant endpoint.

Unlike every other endpoint in this codebase, there is no engine underneath
this one to assert an answer against — what it returns is a language
model's text, which these tests cannot and should not try to grade. What
they prove instead is the HTTP layer around it: auth is required, the
request is bounded, the conversation reaches Ollama in order, and an
Ollama that cannot be reached comes back as a clear 503 rather than a
stack trace or a silent hang.
"""
from unittest.mock import patch

import pytest
from django.urls import reverse

from apps.tutor.ollama import OllamaError


@pytest.mark.django_db
class TestAuthIsRequired:
    def test_signed_out_is_refused(self, api_client):
        r = api_client.post(
            reverse("tutor-ask"), {"question": "What is a sticky end?"}, format="json"
        )
        assert r.status_code == 401

    def test_status_is_also_private(self, api_client):
        r = api_client.get(reverse("tutor-status"))
        assert r.status_code == 401


@pytest.mark.django_db
class TestAvailability:
    def test_enabled_deployment_discloses_ai_limit_before_data_entry(self, auth_client, settings):
        settings.TUTOR_ENABLED = True
        r = auth_client.get(reverse("tutor-status"))
        assert r.status_code == 200
        assert r.data["enabled"] is True
        assert "not a validated scientific result" in r.data["notice"]
        assert r.data["disabled_reason"] == ""

    def test_disabled_deployment_never_calls_ollama(self, auth_client, settings):
        settings.TUTOR_ENABLED = False
        with patch("apps.tutor.views.ask") as mocked:
            r = auth_client.post(
                reverse("tutor-ask"),
                {"question": "What is a sticky end?"},
                format="json",
            )
        assert r.status_code == 503
        assert r.data["code"] == "tutor_disabled"
        mocked.assert_not_called()

    def test_disabled_status_keeps_deterministic_tools_explicitly_available(
        self, auth_client, settings
    ):
        settings.TUTOR_ENABLED = False
        r = auth_client.get(reverse("tutor-status"))
        assert r.status_code == 200
        assert r.data["enabled"] is False
        assert "deterministic" in r.data["disabled_reason"]


@pytest.mark.django_db
class TestRequestIsBounded:
    def test_blank_question_is_rejected(self, auth_client):
        r = auth_client.post(reverse("tutor-ask"), {"question": ""}, format="json")
        assert r.status_code == 400

    def test_question_over_the_limit_is_rejected(self, auth_client):
        r = auth_client.post(reverse("tutor-ask"), {"question": "x" * 2_001}, format="json")
        assert r.status_code == 400

    def test_a_reasonable_question_is_accepted_by_the_serializer(self, auth_client):
        with patch("apps.tutor.views.ask", return_value="ok"):
            r = auth_client.post(
                reverse("tutor-ask"), {"question": "x" * 2_000}, format="json"
            )
        assert r.status_code == 200

    def test_too_much_history_is_rejected(self, auth_client):
        history = [{"role": "user", "content": "hi"}] * 21
        r = auth_client.post(
            reverse("tutor-ask"),
            {"question": "What is a plasmid?", "history": history},
            format="json",
        )
        assert r.status_code == 400

    def test_bad_role_is_rejected(self, auth_client):
        r = auth_client.post(
            reverse("tutor-ask"),
            {"question": "What is a plasmid?",
             "history": [{"role": "system", "content": "hi"}]},
            format="json",
        )
        assert r.status_code == 400

    def test_missing_question_is_rejected(self, auth_client):
        r = auth_client.post(reverse("tutor-ask"), {}, format="json")
        assert r.status_code == 400


@pytest.mark.django_db
class TestTheAnswer:
    def test_returns_what_ollama_said(self, auth_client):
        with patch(
            "apps.tutor.views.ask",
            return_value="A sticky end is a single-stranded overhang.",
        ) as mocked:
            r = auth_client.post(
                reverse("tutor-ask"),
                {"question": "What is a sticky end?", "history": []},
                format="json",
            )
        assert r.status_code == 200
        assert r.data["answer"] == "A sticky end is a single-stranded overhang."
        mocked.assert_called_once()

    def test_history_reaches_ollama_in_order_and_unmodified(self, auth_client):
        history = [
            {"role": "user", "content": "What is codon optimisation?"},
            {"role": "assistant", "content": "Rewriting codon usage for a host."},
        ]
        with patch("apps.tutor.views.ask", return_value="ok") as mocked:
            auth_client.post(
                reverse("tutor-ask"),
                {"question": "Does the protein change?", "history": history},
                format="json",
            )
        called_question, called_history = mocked.call_args[0]
        assert called_question == "Does the protein change?"
        assert called_history == history

    def test_unreachable_ollama_is_a_clear_503_not_a_crash(self, auth_client):
        with patch(
            "apps.tutor.views.ask",
            side_effect=OllamaError("Could not reach Ollama at http://localhost:11434."),
        ):
            r = auth_client.post(
                reverse("tutor-ask"), {"question": "What is a sticky end?"}, format="json"
            )
        assert r.status_code == 503
        assert "Could not reach Ollama" in r.data["detail"]

    def test_real_connection_failure_is_handled(self, auth_client, settings):
        """No mock here: nothing listens on this port, so this exercises the
        real urllib error handling in apps.tutor.ollama rather than trusting
        the mocked path above."""
        settings.OLLAMA_BASE_URL = "http://127.0.0.1:1"
        settings.OLLAMA_TIMEOUT_SECONDS = 2
        r = auth_client.post(
            reverse("tutor-ask"), {"question": "What is a sticky end?"}, format="json"
        )
        assert r.status_code == 503
        assert "Ollama" in r.data["detail"]
