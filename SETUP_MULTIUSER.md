# G-Synth — Multi-user setup

G-Synth now has user accounts. Each visitor signs up, logs in, and gets
a private workspace ("My Projects") that only they can see.

## How it works

| Layer | Technology | Cost |
|---|---|---|
| Hosting | Streamlit Community Cloud | Free |
| Login / accounts | Built in (`utils/auth.py`, PBKDF2-SHA256) | Free |
| Database | SQLAlchemy → SQLite **or** Postgres | Free |

## Two modes

### Demo mode (default — works with zero setup)

If you do nothing, the app runs on a local **SQLite** file. Login and
projects work, **but the data is erased every time Streamlit Cloud
restarts the app** (after a redeploy, or after inactivity). A yellow
"Demo mode" banner reminds you of this.

Fine for trying it out. **Not** fine for real users.

### Persistent mode (recommended — still 100 % free)

Point the app at a free **Supabase Postgres** database. Accounts and
projects then survive forever. No credit card required.

## Connect Supabase (≈ 5 minutes, browser only)

1. Go to **https://supabase.com** → **Start your project** → sign in
   with GitHub.
2. **New project**. Pick a name, set a database password (save it),
   choose the region closest to you, and create. Wait ~2 min for it to
   provision.
3. In the project, open **Project Settings** (gear icon) →
   **Database** → **Connection string** → **URI** tab.
4. Copy the URI. It looks like:
   ```
   postgresql://postgres:[YOUR-PASSWORD]@db.xxxxxxxx.supabase.co:5432/postgres
   ```
   Replace `[YOUR-PASSWORD]` with the password from step 2.
5. Go to **https://share.streamlit.io**, open your G-Synth app →
   **⋮ menu** → **Settings** → **Secrets**.
6. Paste this (one line), then **Save**:
   ```toml
   DATABASE_URL = "postgresql://postgres:YOUR-PASSWORD@db.xxxxxxxx.supabase.co:5432/postgres"
   ```
7. The app reboots automatically. The "Demo mode" banner disappears —
   you're now on persistent storage. The `users` and `projects` tables
   are created automatically on first run.

> Supabase's free tier pauses a project after ~1 week with no activity.
> If that happens, open the Supabase dashboard and click **Resume** —
> your data is intact.

## Running locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Locally it uses a `gsynth.db` SQLite file in the project folder — no
configuration needed. To test against Postgres locally, create
`.streamlit/secrets.toml` with the same `DATABASE_URL` line.

## Files

| File | Role |
|---|---|
| `utils/database.py` | SQLAlchemy models (`User`, `Project`) + CRUD |
| `utils/auth.py` | Password hashing, sign-up / sign-in / sign-out |
| `modules/account.py` | Login gate + account page |
| `modules/my_projects.py` | Per-user saved-projects workspace |

## Security notes

- Passwords are never stored in plain text — only PBKDF2-SHA256 hashes
  (240 000 iterations, per-user 16-byte salt).
- Each database query is scoped to the signed-in `user_id`; users
  cannot read or delete each other's projects.
- Keep `DATABASE_URL` only in Streamlit **Secrets** — never commit it.
  `.streamlit/secrets.toml` is already in `.gitignore`.
