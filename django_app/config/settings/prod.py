"""Production settings — strict, HTTPS-aware, fails fast on missing config.

Every value below is read *without* a default. `django-environ` raises
`ImproperlyConfigured` when the variable is absent, so a misconfigured
deploy dies at startup instead of silently running with a public secret
key, an open host allowlist, or an ephemeral in-container database.
"""
from django.core.exceptions import ImproperlyConfigured

from .base import *  # noqa: F401, F403

DEBUG = False

# ─── Required configuration — no defaults, missing value = hard failure ─────
SECRET_KEY = env("DJANGO_SECRET_KEY")               # noqa: F405
ALLOWED_HOSTS = env.list("ALLOWED_HOSTS")           # noqa: F405
CORS_ALLOWED_ORIGINS = env.list("CORS_ALLOWED_ORIGINS")  # noqa: F405

# Requiring DATABASE_URL is what stops production from falling back to a
# SQLite file on a container's ephemeral disk — where every redeploy would
# silently destroy all accounts and projects.
DATABASES = {"default": env.db_url("DATABASE_URL")}  # noqa: F405

# Belt and braces: the development key is committed to the repository, so
# anyone could forge JWTs with it. Refuse to serve traffic if it leaks into
# a production environment.
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

# ─── HTTPS / hardening ──────────────────────────────────────────────────────
# Behind a TLS-terminating proxy (Render, Fly, Cloudflare)
SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")
SECURE_SSL_REDIRECT = True

# The container's own health probe speaks plain HTTP over loopback and sends
# no X-Forwarded-Proto, so without this exemption SECURE_SSL_REDIRECT answers
# it with a 301 to https://localhost — the probe fails, the orchestrator marks
# the container unhealthy, and it restart-loops forever.
SECURE_REDIRECT_EXEMPT = [r"^api/health/$"]

SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
SECURE_HSTS_SECONDS = 60 * 60 * 24 * 30  # 30 days — bump to 1 year once stable
SECURE_HSTS_INCLUDE_SUBDOMAINS = True
SECURE_HSTS_PRELOAD = False
SECURE_CONTENT_TYPE_NOSNIFF = True
X_FRAME_OPTIONS = "DENY"
