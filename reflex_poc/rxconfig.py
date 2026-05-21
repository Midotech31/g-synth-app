"""Reflex configuration for the G-Synth POC.

`api_url` (the URL the browser uses to reach the FastAPI backend) is only
set explicitly when we *know* it differs from Reflex's own default:

- GitHub Codespaces  → the port-8000 forwarded HTTPS URL
- $API_URL override  → whatever you set
- Local dev          → leave unset, Reflex uses http://localhost:8000
- Reflex Cloud deploy→ leave unset, Reflex Cloud injects the real URL

Forcing `api_url` unconditionally (the previous behaviour) would bake
`localhost:8000` into the production build and break the cloud deploy.
"""
import os

import reflex as rx


def _detect_api_url() -> str | None:
    if override := os.environ.get("API_URL"):
        return override.rstrip("/")
    cs_name = os.environ.get("CODESPACE_NAME")
    if cs_name:
        cs_domain = os.environ.get(
            "GITHUB_CODESPACES_PORT_FORWARDING_DOMAIN", "app.github.dev"
        )
        return f"https://{cs_name}-8000.{cs_domain}"
    return None


_kwargs: dict = {
    "app_name": "gsynth_poc",
    "db_url": "sqlite:///gsynth_poc.db",
    "show_built_with_reflex": False,
}
if (_api_url := _detect_api_url()) is not None:
    _kwargs["api_url"] = _api_url

config = rx.Config(**_kwargs)
