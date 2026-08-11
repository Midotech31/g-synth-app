"""Test settings.

Inherits dev, then lifts the rate limits so ordinary test cases (which hit
register/login repeatedly) don't trip the throttle. The throttling
regression tests re-tighten the rates themselves with `override_settings`,
so the protection is still verified — just not in every unrelated test.
"""
from .dev import *  # noqa: F401, F403

REST_FRAMEWORK = {
    **REST_FRAMEWORK,  # noqa: F405
    "DEFAULT_THROTTLE_RATES": {
        "anon": "10000/hour",
        "user": "10000/hour",
        "register": "10000/hour",
        "login": "10000/min",
        "design": "10000/hour",
        "tutor": "10000/hour",
    },
}

# Fast, deterministic hashing — tests don't need PBKDF2's work factor.
PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]
