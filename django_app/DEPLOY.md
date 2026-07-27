# Deploy the G-Synth API — click-by-click

Goal: a public URL like `https://gsynth-api.onrender.com` that you can open
in a browser, with a database that survives restarts. No terminal, no
Docker, no local Python. Cost: **0 €**.

Two accounts are needed — Supabase (the database) and Render (the server).
Both are free and neither asks for a card.

---

## Part 1 — Database (Supabase) · about 5 minutes

You need this anyway: your **live Streamlit app is currently losing every
account on each restart** because it has no persistent database. The same
Supabase database fixes that app *and* powers this API.

1. Open **https://supabase.com** → **Start your project** → sign in with
   GitHub.
2. **New project**.
   - **Name:** `gsynth`
   - **Database Password:** click *Generate a password*, then **copy it
     somewhere safe** — you need it in step 4 and it is not shown again.
   - **Region:** *Central EU (Frankfurt)* — closest to Algeria.
   - **Create new project**, then wait ~2 minutes while it provisions.
3. Click the green **Connect** button at the top of the page.
4. Choose the **Session pooler** tab and copy the connection string. It
   looks like:

   ```
   postgresql://postgres.abcdefgh:[YOUR-PASSWORD]@aws-0-eu-central-1.pooler.supabase.com:5432/postgres
   ```

   Replace `[YOUR-PASSWORD]` (square brackets included) with the password
   from step 2. **This full line is your `DATABASE_URL`.** Keep it handy.

> Treat it like a password — it grants full access to the database. Never
> paste it into a chat, an issue, or a commit.

### While you're here: fix the live Streamlit app

1. Go to **https://share.streamlit.io** → your G-Synth app → **⋮** →
   **Settings** → **Secrets**.
2. Paste this single line (with your real URL) and **Save**:

   ```toml
   DATABASE_URL = "postgresql://postgres.abcdefgh:YOURPASSWORD@aws-0-eu-central-1.pooler.supabase.com:5432/postgres"
   ```

3. The app restarts by itself. The yellow **"Demo mode"** banner disappears
   — accounts are now permanent.

---

## Part 2 — API server (Render) · about 5 minutes

1. Open **https://render.com** → **Get Started** → sign in with GitHub.
2. In the dashboard: **New +** → **Blueprint**.
3. Pick the repository **`Midotech31/g-synth-app`**. If Render can't see it,
   click *Configure account* and grant access to that repository.
4. Render finds `render.yaml` on its own and shows a service named
   **gsynth-api**. Nothing to configure — except:
5. It asks for **DATABASE_URL**. Paste the Supabase line from Part 1.
6. **Apply** / **Create Blueprint Instance**.

Render now installs dependencies, collects static files, runs the database
migrations, and starts the server. First deploy takes **3–5 minutes**; the
log ends with a line containing `Booting worker`.

Your URL appears at the top of the service page, in the form
**`https://gsynth-api.onrender.com`**.

---

## Check it worked

Open these in your browser:

| URL | What you should see |
|---|---|
| `https://gsynth-api.onrender.com/api/health/` | `{"status": "ok", "service": "gsynth-api"}` |
| `https://gsynth-api.onrender.com/admin/` | The Django login page |

If the health URL returns that JSON, **the API is live**.

### Create your admin account

Service page → **Shell** tab (in the left menu) → type:

```bash
python manage.py createsuperuser
```

Enter your email and a password, then sign in at `/admin/` to browse users
and projects through a web interface.

---

## Things to know about the free tier

- **The service sleeps after 15 minutes of inactivity.** The next visit
  wakes it and takes **~50 seconds** to respond. That's normal on the free
  plan, not a bug. It stays awake while you're using it.
- **Supabase pauses a project after ~1 week with no activity.** If that
  happens, open the Supabase dashboard and click **Resume** — your data is
  intact.
- Both limits disappear on paid plans; neither costs anything to test with.

---

## If something fails

| Symptom | Cause | Fix |
|---|---|---|
| Build fails, log mentions `ImproperlyConfigured: Set the DATABASE_URL environment variable` | DATABASE_URL wasn't entered | Service → **Environment** → add it → **Save**, redeploy |
| Log says `DJANGO_SECRET_KEY is set to the public development key` | The generated key got overwritten | Delete the `DJANGO_SECRET_KEY` variable, redeploy — Render regenerates it |
| `DisallowedHost` errors | The hostname changed | Service → **Environment** → set `ALLOWED_HOSTS` to your exact `xxx.onrender.com` host |
| Health check keeps failing / service restarts in a loop | Database unreachable | Check the Supabase project isn't paused, and that the password in the URL is right |
| `password authentication failed` | `[YOUR-PASSWORD]` was never substituted | Re-copy the connection string and replace the placeholder |

Copy the failing log lines (never the `DATABASE_URL` itself) and I can
diagnose from those.
