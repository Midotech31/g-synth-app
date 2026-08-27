# G-Synth workstation handover

**Prepared:** 2026-08-26 19:30 (Africa/Algiers, UTC+01:00)

**Repository:** `https://github.com/Midotech31/g-synth-app`

**Production branch at audit time:** `main` at `462a20de6b3cc0136f5a37c34b41a37b267fff95`

**Production commit:** `Redesign the scientific workspace (#14)`

**Engine version:** `1.0.0`

**Purpose:** Give a workstation Codex enough product, biological, technical,
deployment, validation, and repository context to continue without access to
the originating chat.

This document contains no credentials. Do not add passwords, personal access
tokens, database URLs, secret keys, or copied environment values to it.

---

## 1. Executive state

G-Synth is a working full-stack laboratory application for one gene-synthesis
and cloning workflow. The current stack is:

```text
React/TypeScript workspace
        ↓ JSON over /api
Django REST API
        ↓ direct function calls
dependency-free Python biology engine
        ↓
PostgreSQL in production / SQLite in development
```

The live services were reachable on 2026-08-26:

- Workspace: `https://gsynth-app.onrender.com/` — HTTP 200.
- API health: `https://gsynth-api-c2p9.onrender.com/api/health/` — returned
  `{"status":"ok","service":"gsynth-api"}`.
- Django admin: `https://gsynth-api-c2p9.onrender.com/admin/`.

Render builds both services from `main`. Supabase provides the persistent
PostgreSQL database through its Session Pooler. Free-tier cold starts are
expected; this deployment does not provide a zero-downtime SLA.

At this handover, every automated functional check passes:

| Layer | Result |
| --- | ---: |
| Focused SSD + PCR compatibility tests | 73 passed |
| Biology engine | 949 passed |
| Django/API | 222 passed |
| Frontend | 47 passed |
| Ruff | passed |
| TypeScript | passed |
| Production frontend build | passed, 502 modules transformed |
| Django migration check | no changes detected |
| Python dependency consistency | no broken requirements |

Non-failing local warnings:

- Pytest could not write its cache in this sandbox.
- Django tests warn that `django_app/staticfiles/` does not exist before
  `collectstatic`; Render creates it during the production build.

Historical note: this baseline had four npm advisories. The 2026-08-27 local
hardening pass resolves them through tested React Router 7 and Vite 8 upgrades;
the current lockfile audits at zero known vulnerabilities. Django was moved
from the unsupported 5.1 line to 5.2.17 after a seven-advisory audit, also with
the complete test cycle.

### Local hardening update — 2026-08-27

The USB working tree now contains an uncommitted product-hardening pass after
the production commit named above. It adds one structured preflight contract
and stable diagnostic codes across optimisation, SSD/assembly, PCR, cloning,
and verification; server-generated provenance persisted on projects and
included in exports; full-duplex validation for supplied inserts; explicit
five-state sequencing verdicts with a coverage map; Guided/Expert modes;
session-restored user-scoped drafts; a cloning bench worksheet; and scientific
limitations, references, worked examples, and a bench-user study protocol in
`docs/`.

Database change: apply `projects.0003_project_provenance` with the normal
`python manage.py migrate` deployment step. Provenance is read-only through the
API; user-editable project JSON is never accepted as an authoritative audit
record.

Validated local state after this pass:

| Layer | Result |
| --- | ---: |
| Biology engine | 1,055 passed |
| Generated/randomized molecular invariants | 70 passed (included above) |
| Engine coverage | 97.21% (95% gate passed) |
| Django/API | 230 passed |
| Frontend | 54 passed |
| Ruff | passed |
| TypeScript | passed |
| Production frontend build | passed, 485 modules transformed |
| Django migration/system checks | no drift; no issues |
| npm audit | 0 known vulnerabilities |
| Python dependency audit | 0 known vulnerabilities |
| Mobile/accessibility engineering audit | contrast, target size and reflow checks passed after fixes |

The frontend was validated in a clean temporary Linux `npm ci` installation.
The dependency folder present on the USB was installed for Windows and cannot
load Rolldown's Linux native binding; use `npm ci` on the target operating
system rather than copying `node_modules` between systems.

The local SQLite database was backed up as
`django_app/db.sqlite3.pre-provenance-20260827.bak`, then migrated successfully
through `projects.0003_project_provenance`. The public Render workspace, SPA
route rewrite, health endpoint, public enzyme catalogue, and unauthenticated
API guard were smoke-tested successfully. The hardened working tree has not
been committed/pushed, so the production PostgreSQL migration is still pending
the authorised deployment; Render's start command will apply it on that deploy.

The dedicated audit is recorded in
`docs/ACCESSIBILITY_MOBILE_AUDIT_2026-08-27.md`. The independent 5–8 scientist
study and physical design-to-sequencing runs remain external evidence gates;
their protocols explicitly remain marked pending and must not be represented
as completed.

---

## 2. Start here on the workstation

Clone the production branch and prove its identity:

```powershell
git clone https://github.com/Midotech31/g-synth-app.git
Set-Location g-synth-app
git switch main
git pull --ff-only
git rev-parse HEAD
```

At the date above, the last command should print:

```text
462a20de6b3cc0136f5a37c34b41a37b267fff95
```

If `main` has advanced, inspect the intervening commits rather than resetting
it to this snapshot.

Create the Python environment on Windows:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
python -m pip install -r django_app\requirements.txt
```

Install the exact frontend dependency lock:

```powershell
Set-Location frontend
npm ci
Set-Location ..
```

Run the baseline before changing code:

```powershell
.\.venv\Scripts\python.exe -m ruff check .
.\.venv\Scripts\python.exe -m pytest gsynth_engine\tests -q
Set-Location django_app
..\.venv\Scripts\python.exe -m pytest -q
..\.venv\Scripts\python.exe manage.py makemigrations --check --dry-run --settings=config.settings.test
Set-Location ..\frontend
npm run typecheck
npm test
npm run build
Set-Location ..
```

To run locally without Docker, use two terminals:

```powershell
# Terminal 1, repository root
Set-Location django_app
..\.venv\Scripts\python.exe manage.py migrate --settings=config.settings.dev
..\.venv\Scripts\python.exe manage.py runserver --settings=config.settings.dev
```

```powershell
# Terminal 2
Set-Location frontend
npm run dev
```

Open `http://localhost:5173`. The frontend proxies `/api` to
`http://127.0.0.1:8000` by default. Docker Compose under `django_app/` also
starts PostgreSQL, Redis, the API, frontend, and a local Ollama service.

---

## 3. What the product is for

G-Synth automates Mohamed Merzoug's laboratory workflow for bacteriocin and
enterocin genes expressed in *E. coli* from pET vectors. It is not intended
to become a generic sequence editor.

The primary flow is:

```text
gene → Optimise → Design → PCR → Clone → Check
                                      ↘ Compare
```

- **Optimise** changes codons for the expression host without changing the
  protein and avoids selected restriction sites.
- **Design** builds the SSD cassette and Merzoug PCR-free oligo assembly.
- **PCR** designs conventional or cloning primers, including enzyme tails,
  digest simulation, primer-quality checks, and reading-frame warnings.
- **Clone** cuts a vector, ligates the insert, re-checks the seams and ORF,
  and draws the plasmid.
- **Check** covers ligation quantities, sequencing primers, text reads, and
  `.ab1` chromatograms with quality-aware difference calls.
- **Compare** performs DNA or protein pairwise alignment.
- **Projects** saves designs/plasmids per user and reopens them.
- **Learn** talks to Ollama when one is configured.

The scientific quality goal is comparable to SnapGene/Geneious while
preserving the lab's custom Merzoug assembly. The differentiator is that
downloads are blocked until the design re-ligates in silico to the intended
construct on both strands.

---

## 4. Legacy executable audit: cloning primers and SSD

### 4.1 Reference and method

The user supplied this reference executable:

```text
Original path: C:\Users\dell\Desktop\G-Synth_2025_4_0.exe
Size: 121,718,276 bytes
SHA-256: 46C56A17DBAF7D99BB6552B19D3631861646847883BFB41F5284C4EBD63ECE9C
File modified: 2025-07-04
```

It is a PyInstaller application built with Python 3.12. The audit extracted
the embedded Python code objects from the supplied executable, identified the
actual primer and SSD functions, reconstructed only those pure calculation
functions, and executed them with the same input molecules as the modern
engine. This is stronger evidence than comparing the current code with the
old `modules/` copy, which is not identical to the executable in one constant.

The GUI itself was not treated as a scientific oracle. The comparison was of
the calculations embedded inside the provided binary.

### 4.2 Overall verdict

The core SSD logic is respected for the validated NdeI/XhoI workflow. Default
and coding SSD designs matched the executable nucleotide for nucleotide in
both forward and reverse oligos.

The modern implementation also corrects real defects in the old generic
enzyme handling and improves PCR thermodynamics and validation. Those changes
should remain.

There is one resolved behavioural compatibility decision, one unresolved
behavioural decision, and one exact DNA difference:

1. **Resolved 2026-08-27:** cloning PCR now defaults to the legacy-safe NdeI
   behaviour: when `CATATG` supplies the start codon, a template-leading
   `ATG` is omitted. The PCR form exposes an explicit opt-in to keep both
   codons and add an N-terminal methionine.
2. Modern SSD accepts a sequence without `ATG` even when `is_coding=True`.
   The executable rejected it.
3. The Factor Xa site encodes the same IEGR peptide but uses different
   synonymous codons.

These are documented decisions, not changes made during this audit.

### 4.3 SSD logic: exact findings

The legacy non-coding cassette order is:

```text
[left enzyme remainder]
[ATG, unless NdeI supplies it]
[GSS left linker]
[6×His]
[SSG right linker]
[optional protease site]
[input sequence]
[right enzyme remainder]
```

The reverse oligo mirrors those pieces in reverse-complement order. The
modern default preserves that order exactly. There is no extra linker after
the input sequence.

| SSD rule | Legacy executable | Modern engine | Verdict |
| --- | --- | --- | --- |
| Default non-coding NdeI/XhoI | Builds the cassette above | Exact same forward and reverse oligos | Preserved |
| Coding NdeI/XhoI | Removes a leading gene `ATG`; NdeI supplies it | Exact same forward and reverse oligos when input begins `ATG` | Preserved |
| Stop removal | Truncates at the first in-frame `TAA`, `TAG`, or `TGA` | Same | Preserved |
| Sequence validation | A/C/G/T only after whitespace cleanup | Same base rule through `validate_dna` | Preserved/improved errors |
| His tag | `CACCACCACCACCACCAC` | Same | Preserved |
| Linkers | `GGTTCTTCT` and `TCTTCTGGT` | Same | Preserved |
| Protease placement | Before the input sequence | Same | Preserved |
| Optional cassette parts | Always present in the old non-coding path | Modern UI/engine can omit His and/or linkers | Intentional extension |
| Coding input without `ATG` | Rejected: “Coding sequence must start with ATG.” | Accepted | Compatibility gap |
| Alternative enzyme ends | Reuses one old remainder pair on both sides; wrong for several enzymes | Derives left/right ends from cut geometry | Correct modern fix |
| Tm | Whole-oligo legacy consensus calculation | SantaLucia nearest-neighbour under annealing conditions | Correct modern upgrade |

The alternative-enzyme difference must not be “fixed” back to the executable.
The executable stored `cut_forward`/`cut_reverse` strings once per enzyme and
reused them at either side of an insert. Left and right molecular roles are
not interchangeable. The modern `left_remainders()` and
`right_remainders()` derive the correct ends from `cut_top`/`cut_bottom` and
are tested for all 109 known enzymes.

Factor Xa exact sequences:

```text
Legacy executable: ATCGAGGGAAGG  → IEGR
Modern engine:     ATCGAAGGTCGT  → IEGR
```

This is peptide-equivalent but not DNA-equivalent. Keep the modern sequence
unless the lab's historical ordered oligos, codon preference, or written SOP
requires the exact legacy DNA. If exact backward reproducibility matters,
make the choice explicit and add a golden test before changing the constant.

### 4.4 Cloning-primer logic: exact findings

The executable's `design_cloning_primers` did the following:

- Selected fixed forward and reverse annealing footprints from the sequence.
- Used a configurable 5' prefix, default `TGCATC`.
- Added the full enzyme recognition site after the prefix.
- For forward NdeI and a gene beginning `ATG`, removed the gene's first three
  bases because `CATATG` already supplies the start codon.
- Built the reverse primer as prefix + right recognition site + reverse
  complement of the chosen final footprint.
- Calculated Tm over the entire tailed oligo.

The modern engine in `gsynth_engine/pcr.py`:

- Chooses 18–30 nt annealing footprints toward a 60 °C target.
- Uses a six-base `GCTAGC` clamp by default so enzymes cut efficiently near
  fragment ends.
- Reports annealing-part Tm separately from whole-oligo Tm and derives Ta from
  the annealing part, which is the physically relevant first-cycle value.
- Defaults to `use_site`: when the left site supplies `ATG`, the forward
  primer anneals at codon two so the expressed protein starts with one
  methionine. `keep_both` remains an explicit opt-in.
- Simulates the PCR product and digest, checks internal sites, verifies sticky
  ends, reports primer-quality problems, and can send the predigested insert
  into Clone.

Matched NdeI/XhoI example:

```text
Legacy forward: TGCATCCATATGAAAGGTGAAGAATTGTT
Modern forward: GCTAGCCATATGAAAGGTGAAGAATTGTTCACCG

Legacy reverse: TGCATCCTCGAGTTTCAGGGTCAGTTTACCGT
Modern reverse: GCTAGCCTCGAGTTTCAGGGTCAGTTTACCGT
```

The length/prefix differences are intentional modernisations. The start-codon
choice is explicit in the PCR interface: the default “Use NdeI's ATG” matches
the executable and avoids Met-Met; “Keep both ATGs” preserves the former
modern behaviour when an extra N-terminal methionine is deliberate. Engine,
digest, reading-frame, API, client, and round-trip tests protect the default.

This interpretation is independently supported by NEB's NdeI definition
(`CA/TATG`) and by a published expression-vector description of the NdeI
`CATATG` site as harbouring the initiation codon:

- https://www.neb.com/en-us/products/r0111-ndei
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4216109/

For SSD, either enforce the documented `is_coding=True` contract (must begin
with `ATG`) or rename/reword the flag so “coding” means only “do not add the
expression cassette”. The former matches the executable and current UI copy.

---

## 5. Architecture and ownership boundaries

```text
gsynth_engine/                 Scientific logic; standard-library runtime only
gsynth_engine/tests/           Molecular claims and golden examples
django_app/apps/design/        Thin HTTP adapters for engine functions
django_app/apps/accounts/      User model, JWT authentication, revocation
django_app/apps/projects/      Per-user saved designs and plasmids
django_app/apps/sequences/     SnapGene/GenBank/FASTA parsing and validation
django_app/apps/tutor/         Optional Ollama adapter
frontend/src/pages/            Workspace screens
frontend/src/components/       Shared scientific/UI renderers
frontend/src/state/            Navigation-lifetime workspace state
tools/                         Generators for committed data/assets
app.py, modules/, utils/       Superseded Streamlit application, retained
```

The most important repository rule is: **no biology in Django views,
serializers, or React components**. Put it in `gsynth_engine`, test the
molecule there