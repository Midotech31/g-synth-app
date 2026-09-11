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


PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]


TUTOR_ENABLED = True
