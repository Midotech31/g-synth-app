from datetime import timedelta

from .base import *  # noqa: F401, F403

DEBUG = True
ALLOWED_HOSTS = ["*"]


CORS_ALLOWED_ORIGINS = [
    "http://localhost:5173",
    "http://localhost:3000",
    "http://127.0.0.1:5173",
]
CORS_ALLOW_ALL_ORIGINS = True


SIMPLE_JWT = {**SIMPLE_JWT, "ACCESS_TOKEN_LIFETIME": timedelta(hours=8)}  # noqa: F405


STORAGES = {
    **STORAGES,  # noqa: F405
    "staticfiles": {"BACKEND": "django.contrib.staticfiles.storage.StaticFilesStorage"},
}


LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {"console": {"class": "logging.StreamHandler"}},
    "loggers": {
        "django": {"handlers": ["console"], "level": "INFO"},
    },
}
