from django.core.exceptions import ImproperlyConfigured

from .base import *  # noqa: F401, F403

DEBUG = False


SECRET_KEY = env("DJANGO_SECRET_KEY")               # noqa: F405
ALLOWED_HOSTS = env.list("ALLOWED_HOSTS")           # noqa: F405


def _https_origin(value):

    origin = value.strip().rstrip("/")
    return origin if "://" in origin else f"https://{origin}"


CORS_ALLOWED_ORIGINS = [
    _https_origin(origin)
    for origin in env.list("CORS_ALLOWED_ORIGINS")  # noqa: F405
]


DATABASES = {"default": env.db_url("DATABASE_URL")}  # noqa: F405


TUTOR_ENABLED = env.bool("TUTOR_ENABLED", default=False)  # noqa: F405


if SECRET_KEY == INSECURE_DEV_SECRET_KEY:           # noqa: F405
    raise ImproperlyConfigured(
        "DJANGO_SECRET_KEY is set to the public development key. "
        "Generate a unique one, e.g. "
        "python -c 'from django.core.management.utils import get_random_secret_key;"
        " print(get_random_secret_key())'"
    )

if "*" in ALLOWED_HOSTS:
    raise ImproperlyConfigured(
        "ALLOWED_HOSTS must not contain '*' in production — list the real hostnames."
    )


SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")
SECURE_SSL_REDIRECT = True


SECURE_REDIRECT_EXEMPT = [r"^api/health/$"]

SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
SECURE_HSTS_SECONDS = 60 * 60 * 24 * 365
SECURE_HSTS_INCLUDE_SUBDOMAINS = True
SECURE_HSTS_PRELOAD = True
SECURE_CONTENT_TYPE_NOSNIFF = True
X_FRAME_OPTIONS = "DENY"
