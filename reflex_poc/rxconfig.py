"""Reflex configuration for the G-Synth POC.

Auto-detects GitHub Codespaces so the same code runs locally and in the
cloud without manual env tweaks.
"""
import os

import reflex as rx


def _detect_api_url() -> str:
    """Return the backend URL the browser should connect to.

    - In GitHub Codespaces, the backend port (8000) is forwarded to a
      public-per-codespace HTTPS URL like
      `https://<name>-8000.app.github.dev`. We assemble that here.
    - Locally, the browser hits the dev server directly on localhost.
    - Override via $API_URL if you need anything custom.
    """
    if override := os.environ.get("API_URL"):
        return override.rstrip("/")
    cs_name = os.environ.get("CODESPACE_NAME")
    cs_domain = os.environ.get(
        "GITHUB_CODESPACES_PORT_FORWARDING_DOMAIN", "app.github.dev"
    )
    if cs_name:
        return f"https://{cs_name}-8000.{cs_domain}"
    return "http://localhost:8000"


config = rx.Config(
    app_name="gsynth_poc",
    db_url="sqlite:///gsynth_poc.db",
    show_built_with_reflex=False,
    api_url=_detect_api_url(),
)
