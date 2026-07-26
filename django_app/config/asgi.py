"""ASGI entry point (uvicorn / daphne).

Loads `.env` before resolving DJANGO_SETTINGS_MODULE — see wsgi.py.
"""
import os
from pathlib import Path

from django.core.asgi import get_asgi_application

_ENV_FILE = Path(__file__).resolve().parent.parent / ".env"
if _ENV_FILE.exists():
    import environ

    environ.Env.read_env(_ENV_FILE)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings.dev")
application = get_asgi_application()
