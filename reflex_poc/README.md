# G-Synth — Reflex POC

Proof-of-concept that replaces Streamlit with **Reflex** (Python → React +
FastAPI). Demonstrates the architecture for a real multi-user G-Synth:

- Email / password authentication (PBKDF2-SHA256)
- Per-user database (SQLite for the POC, swap for Postgres in prod)
- One working bioinformatics module — CRISPR sgRNA designer
- "Save to dashboard" → reopen / delete saved projects
- Interactive UI (no `st.rerun()` page reloads — Reflex compiles to React)

The existing Streamlit app on `main` is **not affected**. This lives in
its own subdirectory so you can compare side-by-side.

## Try it in one click — GitHub Codespaces

You don't need Python, Node, or anything else installed locally.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?repo=Midotech31/g-synth-app&ref=claude/reflex-poc&devcontainer_path=.devcontainer/reflex-poc/devcontainer.json)

1. Click the badge above (signs you in to GitHub if needed).
2. Pick the smallest "2-core / 8 GB" machine and **Create codespace**.
3. Wait ~3 min — the container builds, installs Reflex, runs `reflex init`,
   and starts `reflex run` in the background.
4. When port 3000 finishes forwarding, Codespaces auto-opens the preview
   tab. The very first React compile is silent for ~2 min before the page
   appears — that's normal. Tail `/tmp/reflex.log` in the terminal if you
   want a heartbeat.
5. Sign up with any email + password (≥ 8 chars), the data lives in the
   in-container SQLite for the lifetime of the codespace.

When you're done, stop the codespace from
[github.com/codespaces](https://github.com/codespaces) — free tier
includes 60 hours/month.

## What this is, what it isn't

| | |
|---|---|
| ✅ **Is** | A working end-to-end demo (login → design → save → reopen) you run locally to evaluate Reflex |
| ✅ **Is** | The migration target architecture if you decide to move off Streamlit |
| ❌ **Isn't** | A full port of the 17 G-Synth modules — just the CRISPR designer |
| ❌ **Isn't** | Production-grade (SQLite for the demo, simplified Doench scoring, no email verification, no rate limiting) |

## Run locally

```bash
# from the repo root
cd reflex_poc

# Optional but recommended: virtualenv
python -m venv .venv && source .venv/bin/activate

pip install -r requirements.txt
reflex init --template blank   # one-time: pulls frontend assets
reflex run                      # builds React, starts at http://localhost:3000
```

First run takes 2-3 min (React compilation). Subsequent runs are fast.

## Architecture

```
reflex_poc/
├── rxconfig.py                       Reflex config (db_url, app_name)
├── requirements.txt                  reflex >=0.9
└── gsynth_poc/
    ├── gsynth_poc.py                 entry: rx.App() + routes
    ├── models.py                     SQLModel tables (User, Project)
    ├── auth.py                       AuthState — sign in / sign up / sign out
    ├── crispr_logic.py               pure-Python guide design (no UI deps)
    ├── state.py                      ProjectState — design, save, load
    ├── layout.py                     shared navbar + page wrapper
    └── pages/
        ├── index.py                  /          landing
        ├── login.py                  /login     sign-in form
        ├── signup.py                 /signup    sign-up form
        ├── dashboard.py              /dashboard saved-projects list
        ├── crispr.py                 /crispr    CRISPR designer
        └── project_detail.py         /projects/[id] read-only saved run
```

## Routes

| Route | Auth required | What you see |
|---|---|---|
| `/` | No | Hero + sign-in/up CTAs (or "Go to dashboard" if logged in) |
| `/login` | No | Sign-in form |
| `/signup` | No | Account creation |
| `/dashboard` | Yes | List of saved CRISPR runs |
| `/crispr` | Yes | Sequence input, NGG-PAM scan, save button |
| `/projects/<id>` | Yes (owner only) | Frozen view of a saved run |

## Storage

- SQLite file `gsynth_poc.db` in the working directory.
- Two tables: `user`, `project`.
- Passwords hashed with PBKDF2-SHA256 (200k iterations, 16-byte salt).
- Project payload is JSON in a TEXT column — extends naturally to other
  modules (primers, alignment, etc.) without schema changes.

## Migrating to Supabase / production

Swap two things in `auth.py` and `state.py`:

1. Replace `_hash_password` / `_verify_password` with Supabase Auth SDK
   calls.
2. Replace `with rx.session() as session:` blocks with the Supabase
   client (or keep SQLModel and point `db_url` at the Supabase Postgres
   connection string).

The rest of the app (UI, pages, routing, state) doesn't change.

## Compared to the Streamlit app

| | Streamlit (main) | Reflex POC |
|---|---|---|
| Concurrent users | Session leakage | True isolation |
| Login / auth | None | Native |
| Save a project | No | Yes |
| Reopen later | No | Yes |
| Page reloads on every interaction | Yes | No |
| Custom URL routes | No | Yes |
| Self-host free | Streamlit Cloud | Render / Fly / self |
| Lines of code (this one module) | ~600 (`crispr_designer.py`) | ~150 (`crispr_logic.py` + `pages/crispr.py`) |
