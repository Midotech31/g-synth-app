"""Development settings — verbose, permissive, SQLite by default."""
from .base import *  # noqa: F401, F403

DEBUG = True
ALLOWED_HOSTS = ["*"]

# Every origin is fine in dev — the React dev server usually runs on :5173
CORS_ALLOWED_ORIGINS = [
    "http://localhost:5173",
    "http://localhost:3000",
    "http://127.0.0.1:5173",
]
CORS_ALLOW_ALL_ORIGINS = True

# Short JWT lifetimes are annoying in dev; extend them
from datetime import timedelta  # noqa: E402

SIMPLE_JWT["ACCESS_TOKEN_LIFETIME"] = timedelta(hours=8)   # noqa: F405

# Log SQL queries to the console for debugging
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {"console": {"class": "logging.StreamHandler"}},
    "loggers": {
        "django": {"handlers": ["console"], "level": "INFO"},
    },
}
