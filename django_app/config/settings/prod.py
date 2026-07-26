"""Production settings — strict, HTTPS-aware, DATABASE_URL required."""
from .base import *  # noqa: F401, F403

DEBUG = False

# ALLOWED_HOSTS and CORS_ALLOWED_ORIGINS must come from the environment
# in production — no defaults, so a missing env var fails loudly at startup.
ALLOWED_HOSTS = env("ALLOWED_HOSTS")  # noqa: F405
CORS_ALLOWED_ORIGINS = env("CORS_ALLOWED_ORIGINS")  # noqa: F405

# Behind a TLS-terminating proxy (Render, Fly, Cloudflare)
SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")
SECURE_SSL_REDIRECT = True
SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
SECURE_HSTS_SECONDS = 60 * 60 * 24 * 30  # 30 days — bump to 1 year once stable
SECURE_HSTS_INCLUDE_SUBDOMAINS = True
SECURE_HSTS_PRELOAD = False
SECURE_CONTENT_TYPE_NOSNIFF = True
X_FRAME_OPTIONS = "DENY"
