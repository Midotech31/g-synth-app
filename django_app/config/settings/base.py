"""
Base Django settings for G-Synth web platform.

Environment-driven — read every deployment-specific value through
`environ`. dev.py and prod.py inherit from this file and only override
what actually differs.

IMPORTANT — how defaults work here:
    Defaults are supplied at the *call site* (`env("X", default=...)`),
    never on the `Env()` constructor. A schema-level default silently
    satisfies every later read, which makes it impossible for prod.py to
    demand a value. Keeping defaults at the call site lets prod.py drop
    them and fail loudly on a missing env var.
"""
import sys
from datetime import timedelta
from pathlib import Path

import environ

BASE_DIR = Path(__file__).resolve().parent.parent.parent

# gsynth_engine lives at the repository root, one level above django_app, so
# the design logic is shared with the Streamlit app rather than duplicated.
REPO_ROOT = BASE_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Only DEBUG carries a schema default — it is not security-sensitive and
# every settings module overrides it explicitly anyway.
env = environ.Env(DEBUG=(bool, False))
environ.Env.read_env(BASE_DIR / ".env")

# The value below is public (it is committed to the repository). It exists
# so a fresh checkout runs without setup. prod.py refuses to boot with it.
INSECURE_DEV_SECRET_KEY = "dev-insecure-do-not-use-in-production-32chars-min"

SECRET_KEY = env("DJANGO_SECRET_KEY", default=INSECURE_DEV_SECRET_KEY)
DEBUG = env("DEBUG")
ALLOWED_HOSTS = env.list("ALLOWED_HOSTS", default=["localhost", "127.0.0.1"])

# ─────────────────────────────────────────────────────────────────────────────
# Apps
# ─────────────────────────────────────────────────────────────────────────────
DJANGO_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
]
THIRD_PARTY_APPS = [
    "rest_framework",
    # Persists refresh tokens so they can actually be revoked — required for
    # logout and for cutting off stolen tokens on password change.
    "rest_framework_simplejwt.token_blacklist",
    "corsheaders",
]
LOCAL_APPS = [
    "apps.accounts",
    "apps.projects",
    "apps.sequences",
    "apps.design",
]
INSTALLED_APPS = DJANGO_APPS + THIRD_PARTY_APPS + LOCAL_APPS

# ─────────────────────────────────────────────────────────────────────────────
# Middleware
# ─────────────────────────────────────────────────────────────────────────────
MIDDLEWARE = [
    "corsheaders.middleware.CorsMiddleware",       # must precede CommonMiddleware
    "django.middleware.security.SecurityMiddleware",
    "whitenoise.middleware.WhiteNoiseMiddleware",  # static files in prod
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

ROOT_URLCONF = "config.urls"

TEMPLATES = [{
    "BACKEND": "django.template.backends.django.DjangoTemplates",
    "DIRS": [],
    "APP_DIRS": True,
    "OPTIONS": {"context_processors": [
        "django.template.context_processors.request",
        "django.contrib.auth.context_processors.auth",
        "django.contrib.messages.context_processors.messages",
    ]},
}]

WSGI_APPLICATION = "config.wsgi.application"
ASGI_APPLICATION = "config.asgi.application"

# ─────────────────────────────────────────────────────────────────────────────
# Database — SQLite by default so a fresh checkout runs; prod.py *requires*
# DATABASE_URL so a misconfigured deploy can never silently land on an
# ephemeral SQLite file inside a container.
# ─────────────────────────────────────────────────────────────────────────────
DATABASES = {
    "default": env.db_url(
        "DATABASE_URL",
        default=f"sqlite:///{BASE_DIR / 'db.sqlite3'}",
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
# Auth
# ─────────────────────────────────────────────────────────────────────────────
AUTH_USER_MODEL = "accounts.User"
AUTH_PASSWORD_VALIDATORS = [
    {"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
     "OPTIONS": {"min_length": 8}},
    {"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator"},
    {"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator"},
]

# ─────────────────────────────────────────────────────────────────────────────
# DRF + JWT
# ─────────────────────────────────────────────────────────────────────────────
REST_FRAMEWORK = {
    "DEFAULT_AUTHENTICATION_CLASSES": (
        # Subclass of simplejwt's JWTAuthentication that additionally rejects
        # tokens minted before the user's last credential change.
        "apps.accounts.authentication.VersionedJWTAuthentication",
    ),
    "DEFAULT_PERMISSION_CLASSES": (
        "rest_framework.permissions.IsAuthenticated",
    ),
    "DEFAULT_RENDERER_CLASSES": (
        "rest_framework.renderers.JSONRenderer",
    ),
    "DEFAULT_PAGINATION_CLASS": "rest_framework.pagination.PageNumberPagination",
    "PAGE_SIZE": 25,
    "TEST_REQUEST_DEFAULT_FORMAT": "json",
    # Rate limiting. `register` and `login` are unauthenticated and therefore
    # the two endpoints an attacker can hammer for free — they get their own
    # tighter scopes, applied per-view via ScopedRateThrottle.
    "DEFAULT_THROTTLE_CLASSES": (
        "rest_framework.throttling.AnonRateThrottle",
        "rest_framework.throttling.UserRateThrottle",
    ),
    "DEFAULT_THROTTLE_RATES": {
        "anon": "120/hour",
        "user": "2000/hour",
        "register": "5/hour",
        "login": "10/min",
    },
}

SIMPLE_JWT = {
    "ACCESS_TOKEN_LIFETIME": timedelta(minutes=30),
    "REFRESH_TOKEN_LIFETIME": timedelta(days=14),
    "ROTATE_REFRESH_TOKENS": True,
    # Rotation without blacklisting leaves the superseded refresh token valid,
    # which defeats both logout and password-change revocation.
    "BLACKLIST_AFTER_ROTATION": True,
    "AUTH_HEADER_TYPES": ("Bearer",),
    "TOKEN_OBTAIN_SERIALIZER": "apps.accounts.serializers.VersionedTokenObtainPairSerializer",
}

# ─────────────────────────────────────────────────────────────────────────────
# CORS — the React SPA will call the API from another origin in dev
# ─────────────────────────────────────────────────────────────────────────────
CORS_ALLOWED_ORIGINS = env.list("CORS_ALLOWED_ORIGINS", default=[])
CORS_ALLOW_CREDENTIALS = True

# ─────────────────────────────────────────────────────────────────────────────
# Static files
# ─────────────────────────────────────────────────────────────────────────────
STATIC_URL = "/static/"
STATIC_ROOT = BASE_DIR / "staticfiles"
STORAGES = {
    "default": {"BACKEND": "django.core.files.storage.FileSystemStorage"},
    "staticfiles": {"BACKEND": "whitenoise.storage.CompressedManifestStaticFilesStorage"},
}

# ─────────────────────────────────────────────────────────────────────────────
# Misc
# ─────────────────────────────────────────────────────────────────────────────
LANGUAGE_CODE = "en-us"
TIME_ZONE = "UTC"
USE_I18N = True
USE_TZ = True
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"
