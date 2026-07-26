# G-Synth Django API — Phase 0

Backend for the G-Synth rewrite (see `/gsynth_roadmap` artifact for the full
roadmap). This directory is self-contained: it lives alongside the existing
Streamlit app in the repo but does not depend on it and does not touch it.

## What Phase 0 ships

| Capability | Endpoint | Notes |
|---|---|---|
| Health check | `GET /api/health/` | For Docker HEALTHCHECK, Render, load balancers |
| Sign up | `POST /api/auth/register` | Email + name + password, returns the created user |
| Log in | `POST /api/auth/login` | Email + password, returns access + refresh JWT |
| Refresh | `POST /api/auth/refresh` | Rotate an access token from a refresh token |
| Me | `GET /api/auth/me` | Current user's profile (auth required) |
| Change password | `POST /api/auth/change-password` | Requires current password |
| Projects | `GET/POST /api/projects/` | List / create per-user projects |
| Project detail | `GET/PATCH/DELETE /api/projects/<id>/` | Scoped to the owner |
| Admin | `/admin/` | Django admin for staff users |

Every project query is filtered by `request.user` in the ORM — a user
cannot see or modify another user's rows, and the tests enforce this.

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
pytest                         # runs against SQLite by default
```

Covers: registration validation, login (right + wrong password), JWT
protection, change-password, project CRUD, and per-user isolation
(bob's projects are invisible to alice).

## Layout

```
django_app/
├── config/                 # Django project — settings + URLs + WSGI/ASGI
│   ├── settings/{base,dev,prod}.py
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
├── tests/                  # pytest, conftest with API fixtures
├── manage.py
├── requirements.txt
├── pytest.ini
├── pyproject.toml          # ruff config
├── Dockerfile
├── docker-compose.yml
├── .env.example
└── .github/workflows/ci.yml
```

Same schema (`gsynth_user`, `gsynth_project`) as the Streamlit multi-user
layer already deployed, so pointing both at the same `DATABASE_URL`
gives you a clean migration path.

## Configuration reference

All settings read from environment via `django-environ`. See `.env.example`
for the full list. Required in prod:

| Env var | Notes |
|---|---|
| `DJANGO_SECRET_KEY` | 50+ chars, cryptographically random |
| `DATABASE_URL` | Postgres URL (Supabase, Neon, self-hosted) |
| `ALLOWED_HOSTS` | Comma-separated hostnames the app will accept |
| `CORS_ALLOWED_ORIGINS` | Comma-separated origins the frontend runs from |
| `DJANGO_SETTINGS_MODULE` | `config.settings.prod` |

## What ships next (per roadmap)

- **Phase 1** — React SPA + SeqViz plasmid viewer. Consumes this API.
- **Phase 2** — Port the 17 bioinformatics modules as endpoints under
  `/api/modules/*`, with Celery for the long-running ones.
- **Phase 3** — Interactive sequence + plasmid editor.
- **Phase 4** — Sharing, tabs, collaboration polish.

See the shared roadmap artifact for details.
