"""WSGI entry point (gunicorn / uWSGI).

Loads `.env` before resolving DJANGO_SETTINGS_MODULE so the same
configuration contract holds under gunicorn as under manage.py.
"""
import os
from pathlib import Path

from django.core.wsgi import get_wsgi_application

_ENV_FILE = Path(__file__).resolve().parent.parent / ".env"
if _ENV_FILE.exists():
    import environ

    environ.Env.read_env(_ENV_FILE)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings.dev")
application = get_wsgi_application()
