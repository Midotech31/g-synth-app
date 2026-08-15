"""Regression tests for the security review of commit 907990f.

Each test pins one finding. They are deliberately behavioural — they boot
the real production settings module in a subprocess, or exercise the real
token lifecycle — rather than asserting on the contents of a settings dict,
so they still fail if the protection is removed some other way.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest
from django.core.cache import cache
from django.urls import reverse
from rest_framework_simplejwt.token_blacklist.models import BlacklistedToken

BASE_DIR = Path(__file__).resolve().parent.parent

# base.py reads BASE_DIR/.env, whose values would satisfy the very variables
# these tests deliberately withhold. Skip rather than report a false pass.
DOTENV_PRESENT = (BASE_DIR / ".env").exists()

_BOOT_SNIPPET = """
import django
django.setup()
from django.conf import settings
print("BOOTED")
print("SECRET_PREFIX", settings.SECRET_KEY[:12])
print("ALLOWED_HOSTS", settings.ALLOWED_HOSTS)
print("DB_ENGINE", settings.DATABASES["default"]["ENGINE"])
"""

_COMPLETE_PROD_ENV = {
    "DJANGO_SETTINGS_MODULE": "config.settings.prod",
    "DJANGO_SECRET_KEY": "u" * 60,
    "ALLOWED_HOSTS": "gsynth.example.org",
    "CORS_ALLOWED_ORIGINS": "https://gsynth.example.org",
    "DATABASE_URL": "postgresql://u:p@db.example.org:5432/gsynth",
}


def boot_prod(**overrides) -> subprocess.CompletedProcess:
    """Import config.settings.prod in a clean subprocess.

    `overrides` may set a variable, or remove one by passing None.
    """
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("DJANGO_", "ALLOWED_", "CORS_", "DATABASE_"))}
    env.update(_COMPLETE_PROD_ENV)
    for key, value in overrides.items():
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value
    return subprocess.run(
        [sys.executable, "-c", _BOOT_SNIPPET],
        cwd=BASE_DIR, env=env, capture_output=True, text=True, timeout=90,
    )


@pytest.mark.skipif(DOTENV_PRESENT, reason="a local .env would supply the withheld variables")
class TestProductionRefusesUnsafeConfig:
    """Findings 1-3: prod.py claimed to fail loudly on missing config, but
    schema-level defaults meant it silently ran with a public secret key,
    ALLOWED_HOSTS=['*'] and an ephemeral SQLite database."""

    def test_complete_config_boots(self):
        r = boot_prod()
        assert "BOOTED" in r.stdout, r.stderr
        assert "postgresql" in r.stdout

    def test_missing_secret_key_is_fatal(self):
        r = boot_prod(DJANGO_SECRET_KEY=None)
        assert r.returncode != 0
        assert "DJANGO_SECRET_KEY" in (r.stderr + r.stdout)

    def test_public_dev_secret_key_is_rejected(self):
        from config.settings.base import INSECURE_DEV_SECRET_KEY
        r = boot_prod(DJANGO_SECRET_KEY=INSECURE_DEV_SECRET_KEY)
        assert r.returncode != 0
        assert "development key" in (r.stderr + r.stdout)

    def test_missing_allowed_hosts_is_fatal(self):
        r = boot_prod(ALLOWED_HOSTS=None)
        assert r.returncode != 0
        assert "ALLOWED_HOSTS" in (r.stderr + r.stdout)

    def test_wildcard_allowed_hosts_is_rejected(self):
        r = boot_prod(ALLOWED_HOSTS="*")
        assert r.returncode != 0
        assert "ALLOWED_HOSTS" in (r.stderr + r.stdout)

    def test_missing_database_url_is_fatal(self):
        """Without this, prod silently falls back to SQLite on a container's
        ephemeral disk and loses every account on redeploy."""
        r = boot_prod(DATABASE_URL=None)
        assert r.returncode != 0, (
            "prod booted without DATABASE_URL — it would silently use SQLite:\n"
            + r.stdout
        )
        assert "sqlite" not in r.stdout.lower()

    def test_missing_cors_origins_is_fatal(self):
        r = boot_prod(CORS_ALLOWED_ORIGINS=None)
        assert r.returncode != 0

    def test_render_host_is_normalised_to_https_origin(self):
        snippet = _BOOT_SNIPPET + "\nprint('CORS', settings.CORS_ALLOWED_ORIGINS)\n"
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("DJANGO_", "ALLOWED_", "CORS_", "DATABASE_"))}
        env.update(_COMPLETE_PROD_ENV)
        env["CORS_ALLOWED_ORIGINS"] = "gsynth-app.onrender.com"
        r = subprocess.run(
            [sys.executable, "-c", snippet], cwd=BASE_DIR, env=env,
            capture_output=True, text=True, timeout=90,
        )
        assert r.returncode == 0, r.stderr
        assert "https://gsynth-app.onrender.com" in r.stdout


@pytest.mark.skipif(DOTENV_PRESENT, reason="a local .env would supply the withheld variables")
def test_health_endpoint_is_exempt_from_ssl_redirect():
    """Finding 6: SECURE_SSL_REDIRECT answered the container's own plain-HTTP
    health probe with a 301 to https://localhost, so Docker/Render marked the
    container unhealthy and restart-looped it."""
    snippet = """
import django
django.setup()
from django.test import Client
r = Client().get("/api/health/", HTTP_HOST="gsynth.example.org")
print("HEALTH_STATUS", r.status_code)
r2 = Client().get("/admin/login/", HTTP_HOST="gsynth.example.org")
print("OTHER_STATUS", r2.status_code)
"""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("DJANGO_", "ALLOWED_", "CORS_", "DATABASE_"))}
    env.update(_COMPLETE_PROD_ENV)
    env["DATABASE_URL"] = "sqlite:///health-probe.sqlite3"  # never queried
    r = subprocess.run([sys.executable, "-c", snippet], cwd=BASE_DIR, env=env,
                       capture_output=True, text=True, timeout=90)
    assert "HEALTH_STATUS 200" in r.stdout, f"health probe not exempt:\n{r.stdout}{r.stderr}"
    # Everything else must still be redirected to HTTPS.
    assert "OTHER_STATUS 301" in r.stdout, f"SSL redirect disabled entirely:\n{r.stdout}"


@pytest.mark.django_db
class TestTokenRevocation:
    """Finding 4: changing the password left both the old access token and
    the old refresh token fully usable — the attacker kept access for the
    refresh lifetime (14 days)."""

    def _login(self, api_client, email="v@example.com", password="GoodPass123"):
        api_client.post(reverse("auth-register"), {
            "email": email, "name": "V",
            "password": password, "password2": password,
        })
        r = api_client.post(reverse("auth-login"), {"email": email, "password": password})
        assert r.status_code == 200, r.data
        return r.data["access"], r.data["refresh"]

    def test_access_token_dies_on_password_change(self, api_client):
        access, _refresh = self._login(api_client)
        api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {access}")
        assert api_client.get(reverse("auth-me")).status_code == 200

        r = api_client.post(reverse("auth-change-password"), {
            "current_password": "GoodPass123",
            "new_password": "Rotated456xyz", "new_password2": "Rotated456xyz",
        })
        assert r.status_code == 200

        # Same token, same request — must now be refused.
        assert api_client.get(reverse("auth-me")).status_code == 401

    def test_refresh_token_dies_on_password_change(self, api_client):
        _access, refresh = self._login(api_client, email="r@example.com")
        api_client.credentials(
            HTTP_AUTHORIZATION=f"Bearer {_access}")
        api_client.post(reverse("auth-change-password"), {
            "current_password": "GoodPass123",
            "new_password": "Rotated456xyz", "new_password2": "Rotated456xyz",
        })
        api_client.credentials()
        r = api_client.post(reverse("auth-refresh"), {"refresh": refresh})
        assert r.status_code == 401, "revoked refresh token still mints access tokens"

    def test_new_password_lets_you_back_in(self, api_client):
        access, _ = self._login(api_client, email="back@example.com")
        api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {access}")
        api_client.post(reverse("auth-change-password"), {
            "current_password": "GoodPass123",
            "new_password": "Rotated456xyz", "new_password2": "Rotated456xyz",
        })
        api_client.credentials()
        r = api_client.post(reverse("auth-login"),
                            {"email": "back@example.com", "password": "Rotated456xyz"})
        assert r.status_code == 200
        api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {r.data['access']}")
        assert api_client.get(reverse("auth-me")).status_code == 200

    def test_logout_blacklists_that_refresh_token(self, api_client):
        access, refresh = self._login(api_client, email="lo@example.com")
        api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {access}")
        assert api_client.post(reverse("auth-logout"), {"refresh": refresh}).status_code == 205
        assert BlacklistedToken.objects.count() == 1
        api_client.credentials()
        assert api_client.post(reverse("auth-refresh"), {"refresh": refresh}).status_code == 401

    def test_logout_all_revokes_every_session(self, api_client):
        access_a, refresh_a = self._login(api_client, email="multi@example.com")
        # A second device signs in with the same credentials.
        r = api_client.post(reverse("auth-login"),
                            {"email": "multi@example.com", "password": "GoodPass123"})
        access_b, refresh_b = r.data["access"], r.data["refresh"]

        api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {access_b}")
        assert api_client.post(reverse("auth-logout-all")).status_code == 200

        for token in (access_a, access_b):
            api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {token}")
            assert api_client.get(reverse("auth-me")).status_code == 401
        api_client.credentials()
        for token in (refresh_a, refresh_b):
            assert api_client.post(reverse("auth-refresh"), {"refresh": token}).status_code == 401

    def test_token_without_version_claim_is_refused(self, api_client, user):
        """Fail closed on tokens minted before versioning existed — and on
        any future call site that bypasses VersionedRefreshToken."""
        from rest_framework_simplejwt.tokens import RefreshToken
        token = RefreshToken.for_user(user).access_token   # deliberately unversioned
        assert "ver" not in token.payload
        api_client.credentials(HTTP_AUTHORIZATION=f"Bearer {token}")
        assert api_client.get(reverse("auth-me")).status_code == 401


@pytest.mark.django_db
class TestThrottling:
    """Finding 7: no throttling at all — 10/10 accounts created back to back,
    and unlimited password guesses against login."""

    def _tighten(self, monkeypatch, **rates):
        """Lower the limits for the duration of one test.

        `SimpleRateThrottle.THROTTLE_RATES` is bound to the settings dict at
        class-definition time, so `override_settings` alone would not reach
        it — patch the class attribute directly.
        """
        from rest_framework.throttling import SimpleRateThrottle
        cache.clear()
        for scope, rate in rates.items():
            monkeypatch.setitem(SimpleRateThrottle.THROTTLE_RATES, scope, rate)

    def test_register_is_throttled(self, api_client, monkeypatch):
        self._tighten(monkeypatch, register="2/hour")
        codes = [
            api_client.post(reverse("auth-register"), {
                "email": f"spam{i}@example.com", "name": "S",
                "password": "GoodPass123", "password2": "GoodPass123",
            }).status_code
            for i in range(5)
        ]
        cache.clear()
        assert 429 in codes, f"registration never throttled: {codes}"

    def test_the_expensive_design_endpoints_are_throttled(
        self, auth_client, monkeypatch,
    ):
        """Optimisation, alignment and verification each cost real CPU.

        Left on the default user rate, one signed-in person can ask for more
        compute than the worker has — which on a single-process deployment
        means the app stops answering anyone, without anything malicious
        happening.
        """
        self._tighten(monkeypatch, design="3/hour")
        gene = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGG"
        codes = [
            auth_client.post(reverse("design-optimise"), {"sequence": gene}).status_code
            for _ in range(6)
        ]
        assert 429 in codes, f"never throttled: {codes}"

    def test_the_scope_is_shared_across_the_design_endpoints(
        self, auth_client, monkeypatch,
    ):
        """One budget for all of them, or the limit is trivially sidestepped
        by spreading the load over four routes."""
        self._tighten(monkeypatch, design="3/hour")
        gene = "ATGACAACAAGTAAATTAGGGAAAGGTTTAGGG"

        for _ in range(3):
            auth_client.post(reverse("design-optimise"), {"sequence": gene})
        spent = auth_client.post(reverse("design-align"), {
            "first": gene, "second": gene,
        })
        assert spent.status_code == 429

    def test_reference_lookups_are_not_throttled_with_them(self, api_client):
        """The enzyme and vector tables are lookups, not work."""
        for _ in range(6):
            assert api_client.get(reverse("design-enzymes")).status_code == 200
            assert api_client.get(reverse("design-vectors")).status_code == 200

    def test_login_is_throttled(self, api_client, monkeypatch, user):
        self._tighten(monkeypatch, login="3/min")
        codes = [
            api_client.post(reverse("auth-login"),
                            {"email": user.email, "password": "wrong-guess"}).status_code
            for _ in range(6)
        ]
        cache.clear()
        assert 429 in codes, f"password guessing never throttled: {codes}"
