# G-Synth toward Geneious-class — technical roadmap

*Technical proposal · v0.1 · July 2026*

**Prepared for** Dr. Mohamed Merzoug — Higher School of Biological Sciences of Oran (ESSBO) · Genomics Technology Platform.

📄 **Illustrated version (shareable URL):** https://claude.ai/code/artifact/726ffb26-adb5-407a-902d-cf9b27da0eb8

---

## In one paragraph

The current G-Synth runs seventeen validated bioinformatics modules on top of Streamlit. The science is sound; the framework is the limit. Streamlit rerenders the whole page on every interaction, cannot render a real interactive plasmid map, and gives no first-class way to hold a multi-user project graph. This proposal keeps every line of validated science, and replaces the shell — **Django REST on the server, React with SeqViz on the client, PostgreSQL and Celery underneath**. Five phases, sequenced so that each one ships something usable.

| | |
|---|---|
| Roadmap length | **5 phases** |
| Solo effort | **6 – 9 months** |
| With one frontend collaborator | **3 – 4 months** |
| Hosting to ≈100 users | **0 $ / mo** (Render free + Supabase free + Cloudflare R2) |
| Existing modules ported | **17 / 17** |

## Where we are, where we need to go

| Capability | G-Synth today | SnapGene | Geneious Prime | Target |
|---|---|---|---|---|
| Multi-user accounts | Just shipped | Yes | Yes | Keep + roles |
| Interactive plasmid map | Static image only | Full editor | Full editor | SeqViz + edit |
| Sequence editing (Word-style) | No | Yes | Yes | Phase 3 |
| Multi-tab project workspace | One view at a time | Yes | Yes | Phase 4 |
| Auto-annotation of features | Basic list | Yes, curated DB | Yes | plannotate + curated |
| SnapGene `.dna` read/write | Read only | Native | Read | Read + best-effort write |
| Chromatograms (`.ab1`) | Not supported | Yes | Yes | Read + viewer |
| CRISPR guide design | Doench 2016 | Basic | Plugins | Ported + Azimuth |
| Collaboration & sharing | None | Cloud add-on | Full | Links + read-only view |
| Cost to end user | Free | ≈ 800 $ / yr | ≈ 1 500 $ / yr | Free (public tier) |

## Target architecture

Five moving parts. The browser holds the interactive workspace; the Django service owns identity, project data, and every bioinformatics endpoint; long-running jobs (alignment, docking, feature scanning) run out-of-band on Celery so the request stays snappy; PostgreSQL and object storage sit at the bottom. Everything speaks JSON over HTTPS.

```
┌────────────────────────────────────────────────────────────┐
│                    User's browser                          │
│  React 18 SPA · TypeScript · Vite · Tailwind · shadcn/ui   │
│  SeqViz plasmid viewer · Monaco sequence editor            │
└────────────────────┬───────────────────────────────────────┘
                     │ HTTPS · JWT
┌────────────────────┴───────────────────────────────────────┐
│                Django REST API (Django 5)                  │
│         DRF · dj-rest-auth · biopython · primer3-py        │
└──────┬──────────────────────┬──────────────────────────────┘
       │                      │
┌──────┴─────┐        ┌───────┴────────┐       ┌───────────┐
│  Celery +  │        │   PostgreSQL   │       │  S3 / R2  │
│  Redis     │        │   (Supabase)   │       │  blobs    │
└────────────┘        └────────────────┘       └───────────┘
```

## Tech stack

**Backend** — Django 5, Django REST Framework + JWT (simplejwt), Celery + Redis (long jobs), biopython, snapgene-reader, reportlab / weasyprint (PDF exports).

**Frontend** — React 18 + Vite + TypeScript, Tailwind + shadcn/ui, SeqViz (linear + circular plasmid viewer), TanStack Query (caching/mutations), Zustand (editor/tab state), CodeMirror 6 (sequence editor).

**Data** — PostgreSQL (Supabase free tier), Redis (Celery broker + view caching), Cloudflare R2 (S3-compatible object storage, no egress fees).

**Infrastructure** — Docker + compose, Render (or Fly.io), Cloudflare (DNS + TLS + CDN + Turnstile), Sentry free tier, GitHub Actions CI.

## Roadmap in five phases

### Phase 0 — Foundation (2 weeks)

Scaffold the project, pick a home for it. Nothing science yet — just the spine everything else stands on.

- **Deliverables:** Django 5 project · custom User model · JWT endpoints (login/refresh/register/me/change-password) · Postgres via Supabase · Dockerfile + docker-compose · GitHub Actions running ruff + pytest.
- **Key libraries:** `django, djangorestframework, djangorestframework-simplejwt, django-cors-headers, django-environ, psycopg[binary]`
- **Exit:** A deployed staging URL where `POST /api/auth/register` creates a user in Supabase Postgres and `POST /api/auth/login` returns a JWT.

### Phase 1 — Sequence viewer MVP (4 weeks)

The smallest slice that already feels like a real tool: user logs in, uploads a FASTA or GenBank, sees a linear + circular map, saves it as a project.

- **Deliverables:** React SPA (Vite/TS) · sign-up / sign-in screens · file-upload endpoint (FASTA · GenBank · `.dna`) · SeqViz linear & circular render · project CRUD API + UI · Tailwind design tokens locked.
- **Key libraries:** `react, vite, typescript, tailwindcss, @radix-ui, seqviz, @tanstack/react-query, zustand, react-router`
- **Exit:** Upload the pUC19 GenBank file → see the annotated circular map in the browser → save as project → close tab → sign in from another device → open the same project.

### Phase 2 — Port the seventeen modules (10 weeks)

The biggest phase by volume, the least risky by nature. Each Streamlit module becomes a Django endpoint + a React page, in this order — quick wins first, heavy compute last.

- Batch A (small, stateless): reverse complement · translation · codon optimization · hybridization · ligation calculator · ligation check.
- Batch B (compute + tables): primer generator · CRISPR sgRNA designer · alignment tools · extended synthesis · SSD.
- Batch C (heavy / long-running): plasmid visualizer · functional prediction · in-silico docking.
- **Async pattern:** anything > 5 s of compute is a Celery task; API returns a job id, the SPA polls for completion.
- **Exit:** Every workflow a G-Synth user does today can be done in the new app, with results persisted to their account.

### Phase 3 — Interactive editor — the SnapGene moment (8 weeks)

This is where the tool starts to feel like SnapGene. Click-and-drag on the plasmid map to move a feature. Right-click a selection to annotate it. Type into the sequence, watch the Tm and GC recompute live. This is the hardest phase and the one that justifies the whole rewrite.

- **Deliverables:** sequence editor with keyboard selection, cut/copy/paste, undo/redo · in-place feature annotation · drag-to-resize on the plasmid map · live GC/Tm readouts · optimistic UI with server-side reconciliation.
- **Key libraries:** `codemirror 6, immer, yjs` (optional, sets up Phase 4 collab).
- **Exit:** A user can build a pUC19 → insert cloning entirely in the browser.

### Phase 4 — Collaboration + polish (8 weeks)

Multi-tab session (open several projects at once). Share a project by link — read-only or read-write for a collaborator. Comments pinned to a feature. Export polish (PDF reports, SVG plasmid maps). Performance tuning for large sequences (> 100 kb).

- **Exit:** A PI can share a plasmid design with a student, both edit it from different browsers, and export a clean PDF report of the cloning steps for their lab notebook.

## Real risks

| Severity | Risk | Mitigation |
|---|---|---|
| **High** | Phase 3 frontend depth (editor + interactive map) | Plan Phase 3 with a hired React collaborator (2 months) or split 3.a "read-only interactive" and 3.b "editing". |
| **Medium** | Feature parity with SnapGene is a moving target | Pin the target to the ESSBO / academic workflow — build SnapGene for that workflow, not for everyone. |
| **Medium** | SnapGene `.dna` write is only partially reverse-engineered | Read `.dna`; export GenBank (SnapGene opens it natively). Skip `.dna` write until a user asks. |
| **Medium** | Large plasmids (> 100 kb) render performance | Viewport culling (CodeMirror does this) + level-of-detail on the circular map. Phase 4 concern. |
| **Low** | Hosting cost creep past ~100 users | Instrument from day one (Sentry + PostHog free tiers). Cost alert when leaving free tier. |
| **Low** | Recruiting a React + bioinformatics dev | Interactive editor is 90 % React, 10 % biology. Hire strong React, provide the biology inline. |

## What to do this week

1. **Lock the deployment target.** Render (simple, generous free tier, git-based deploys) is the default. Alternative: Fly.io.
2. **Provision Supabase now.** Free Postgres, connect it to the current Streamlit for continuity (`DATABASE_URL` secret). Same DB becomes Phase 0's Django DB — zero data migration later.
3. **Create the new backend repo — or use `django_app/` in this repo.** Phase 0 scaffolding is in `django_app/`; see its README for the one-command boot.
4. **Green-light Phase 0.** Two focused weeks. If a collaborator will join, that decision belongs at the boundary between Phase 0 and Phase 1.

---

*Revision 0.1 · Assumes a full-time engineer at focused pace. Halve the speed for part-time; double it for a distracted week.*
