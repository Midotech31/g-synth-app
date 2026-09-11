#!/usr/bin/env python

import os
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def load_dotenv() -> None:

    env_file = BASE_DIR / ".env"
    if not env_file.exists():
        return
    try:
        import environ
    except ImportError:
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
