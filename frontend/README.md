# G-Synth frontend — Phase 1

React SPA for the G-Synth workspace. Talks to the Django API in
`../django_app`.

## What Phase 1 does

- **Sign up / sign in** against the JWT API, with automatic access-token
  refresh (concurrent 401s share one refresh, so a page issuing several
  requests can't invalidate its own rotating refresh token).
- **Import a sequence** — drag a FASTA or GenBank file onto the dashboard,
  or use the file picker. The server parses it with biopython and returns
  annotations already shaped for the viewer.
- **Plasmid map** — circular, linear, or both, with features drawn in
  direction-aware arrows, coloured by feature type, and cross-listed in a
  sidebar with coordinates.
- **Projects persist** to your account. Sign out, sign in on another
  machine, the map is still there.

## Run it

Two processes. Backend first:

```bash
cd ../django_app
pip install -r requirements.txt
python manage.py migrate
python manage.py runserver          # http://127.0.0.1:8000
```

Then the frontend:

```bash
npm install
npm run dev                          # http://localhost:5173
```

Open **http://localhost:5173**. Vite proxies `/api` to Django, so the
browser sees a single origin and there is no CORS preflight in
development. Point it elsewhere with `VITE_API_TARGET`:

```bash
VITE_API_TARGET=https://gsynth-api.onrender.com npm run dev
```

## Build

```bash
npm run typecheck    # tsc --noEmit
npm run build        # → dist/
```

## Layout

```
frontend/
├── index.html
├── vite.config.ts          # dev server + /api proxy
└── src/
    ├── main.tsx
    ├── App.tsx             # routes + app shell + route guards
    ├── styles.css          # design tokens and all component styling
    ├── api/client.ts       # fetch wrapper, JWT storage, refresh-on-401
    ├── auth/AuthContext.tsx
    └── pages/
        ├── Login.tsx
        ├── Signup.tsx
        ├── Dashboard.tsx   # project list + drag-and-drop import
        └── Viewer.tsx      # SeqViz map, stats, feature list
```

## Design notes

A deliberate single-theme design — navy chrome around a warm paper
canvas — rather than a light/dark pair, matching the desktop tools this
replaces. Colour carries meaning: feature types have fixed colours
(assigned server-side in `apps/sequences/parsing.py`, so the map and the
feature list can never disagree), and the single teal accent marks the
primary action on each screen.

The roadmap named Tailwind + shadcn/ui. This ships hand-written CSS with
tokens instead: at Phase 1's size it is fewer moving parts, and the tokens
in `styles.css` map 1:1 onto Tailwind theme values if we adopt it in
Phase 3, when the editor UI grows enough to justify a component library.

## Not in Phase 1

Editing. The map is read-only — clicking a feature highlights it, but you
cannot drag, annotate, or type into the sequence yet. That is Phase 3, and
it is the hard one.
