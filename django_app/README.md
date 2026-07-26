# G-Synth Django API — Phase 0

Backend for the G-Synth rewrite (see `/gsynth_roadmap` artifact for the full
roadmap). This directory is self-contained: it lives alongside the existing
Streamlit app in the repo but does not depend on it and does not touch it.

## What Phase 0 ships

| Capability | Endpoint | Notes |
|---|---|---|
| Health check | `GET /api/health/` | For Docker HEALTHCHECK, Render, load balancers |
| Sign up | `POST /api/auth/register` | Email + name + password. Rate-limited (5/hour/IP) |
| Log in | `POST /api/auth/login` | Email + password → access + refresh JWT. Rate-limited (10/min/IP) |
| Refresh | `POST /api/auth/refresh` | Rotate an access token from a refresh token |
| Log out | `POST /api/auth/logout` | Blacklists the supplied refresh token (this device) |
| Log out everywhere | `POST /api/auth/logout-all` | Revokes every session for the account |
| Me | `GET /api/auth/me` | Current user's profile (auth required) |
| Change password | `POST /api/auth/change-password` | Requires current password; **revokes all existing tokens** |
| Projects | `GET/POST /api/projects/` | List / create per-user projects |
| Project detail | `GET/PATCH/DELETE /api/projects/<id>/` | Scoped to the owner |
| Admin | `/admin/` | Django admin for staff users |

Every project query is filtered by `request.user` in the ORM — a user
cannot see or modify another user's rows, and the tests enforce this.

### Token revocation

JWTs are self-contained and normally stay valid until they expire, which
means a stolen token survives a password change — precisely when it must
not. Two mechanisms close that:

- Every token carries a `ver` claim (the user's `token_version`). Changing
  the password bumps the counter, so tokens issued earlier are refused by
  `VersionedJWTAuthentication` on their next request.
- Outstanding refresh tokens are blacklisted, so they can no longer be
  exchanged for fresh access tokens.

Mint tokens through `apps.accounts.tokens.VersionedRefreshToken` — a bare
`RefreshToken.for_user()` omits the claim and the token will be rejected.

## Local dev — one-command boot

```bash
cd django_app
cp .env.example .env

# Option A — Docker (Postgres + Django in one shot)
docker compose up --build
# API at http://localhost:8000/api/health/

# Option B — Native (SQLite, no Postgres needed)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python manage.py migrate
python manage.py createsuperuser
python manage.py runserver
```

## Try it (from another terminal)

```bash
# Register
curl -X POST http://localhost:8000/api/auth/register/ \
    -H "Content-Type: application/json" \
    -d '{"email":"you@example.com","name":"You","password":"GoodPass123","password2":"GoodPass123"}'

# Log in — grab the access token
ACCESS=$(curl -s -X POST http://localhost:8000/api/auth/login/ \
    -H "Content-Type: application/json" \
    -d '{"email":"you@example.com","password":"GoodPass123"}' \
    | python -c "import sys,json;print(json.load(sys.stdin)['access'])")

# Create a project
curl -X POST http://localhost:8000/api/projects/ \
    -H "Content-Type: application/json" -H "Authorization: Bearer $ACCESS" \
    -d '{"name":"Insulin v1","module":"crispr_designer","sequence":"ATGAAACGT"}'
```

## Tests

```bash
pytest                         # 35 tests, runs against SQLite
```

Covers: registration validation, login (right + wrong password), JWT
protection, change-password, project CRUD, and per-user isolation
(bob's projects are invisible to alice).

`tests/test_security.py` holds regression tests for the findings of the
2026-07 review — they boot the real production settings in a subprocess and
exercise the real token lifecycle, so they fail if a protection is removed:
production refusing unsafe config, token revocation on password change and
logout, the health-check SSL exemption, and rate limiting.

## Layout

```
django_app/
├── config/                 # Django project — settings + URLs + WSGI/ASGI
│   ├── settings/{base,dev,prod,test}.py
│   ├── urls.py             # top-level routes
│   └── wsgi.py / asgi.py
├── apps/
│   ├── accounts/           # custom User (email-login) + auth endpoints
│   │   ├── models.py
│   │   ├── serializers.py
│   │   ├── views.py
│   │   └── urls.py
│   └── projects/           # per-user Project model + CRUD viewset
│       ├── models.py
│       ├── serializers.py
│       ├── views.py
│       └── urls.py
├── tests/                  # pytest — API fixtures + security regressions
├── manage.py
├── requirements.txt
├── pytest.ini              # runs against config.settings.test
├── pyproject.toml          # ruff config
├── Dockerfile
├── docker-compose.yml
└── .env.example

../.github/workflows/django-ci.yml    # CI lives at the REPO ROOT — GitHub
                                      # ignores workflow files anywhere else
```

Same schema (`gsynth_user`, `gsynth_project`) as the Streamlit multi-user
layer already deployed, so pointing both at the same `DATABASE_URL`
gives you a clean migration path.

## Configuration reference

All settings read from environment via `django-environ`. `.env` is loaded
before the settings module is chosen, so `DJANGO_SETTINGS_MODULE` works
there too. See `.env.example` for the full list.

**`config.settings.prod` refuses to start** if any of these is missing —
deliberately, so a misconfigured deploy fails at boot instead of silently
running with a public secret key or an ephemeral database:

| Env var | Notes |
|---|---|
| `DJANGO_SECRET_KEY` | 50+ chars, random. Booting with the committed dev key is rejected |
| `DATABASE_URL` | Postgres URL (Supabase, Neon, self-hosted). Without it, prod would fall back to container-local SQLite and lose all data on redeploy |
| `ALLOWED_HOSTS` | Comma-separated hostnames. `*` is rejected in prod |
| `CORS_ALLOWED_ORIGINS` | Comma-separated origins the frontend runs from |
| `DJANGO_SETTINGS_MODULE` | `config.settings.prod` |

Generate a key with:

```bash
python -c "from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())"
```

## What ships next (per roadmap)

- **Phase 1** — React SPA + SeqViz plasmid viewer. Consumes this API.
- **Phase 2** — Port the 17 bioinformatics modules as endpoints under
  `/api/modules/*`, with Celery for the long-running ones.
- **Phase 3** — Interactive sequence + plasmid editor.
- **Phase 4** — Sharing, tabs, collaboration polish.

See the shared roadmap artifact for details.
