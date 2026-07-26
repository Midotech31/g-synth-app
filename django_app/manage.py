#!/usr/bin/env python
"""Django management CLI.

`.env` is loaded *before* DJANGO_SETTINGS_MODULE is resolved, so setting
that variable in `.env` genuinely selects the settings module. Reading it
only inside settings/base.py (as this file used to) would be too late —
the module had already been chosen and the line silently did nothing.
"""
import os
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def load_dotenv() -> None:
    """Populate os.environ from .env without overriding the real environment."""
    env_file = BASE_DIR / ".env"
    if not env_file.exists():
        return
    try:
        import environ
    except ImportError:  # dependencies not installed yet
        return
    environ.Env.read_env(env_file)


def main() -> None:
    load_dotenv()
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "config.settings.dev")
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Django is not installed. Activate the venv and `pip install -r requirements.txt`."
        ) from exc
    execute_from_command_line(sys.argv)


if __name__ == "__main__":
    main()
