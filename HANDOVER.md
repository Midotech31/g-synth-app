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
molecule there, then expose and render the returned result.

When adding an endpoint:

1. Implement and test the engine function.
2. Bound every serializer input.
3. Translate `SequenceError` into a useful 400 response.
4. Add design throttling to algorithmic work.
5. Add an HTTP test that asserts the response equals the engine result.
6. Render it without re-deriving sequences in TypeScript.

The engine is dependency-free at runtime. Biopython is used in tests,
generated-data tools, and Django file parsing, not by the engine's public
calculation path.

---

## 6. Scientific invariants

Treat a failure here as a wrong molecule, not a normal software defect:

1. Re-ligating a Merzoug plan reproduces the intended construct base for
   base on both strands before a download is offered.
2. Terminal sticky ends are measured from the assembled duplex, not copied
   from labels, and match the selected enzymes in sequence and polarity.
3. Codon optimisation never changes the translated protein.
4. Recutting a recombinant plasmid with the same pair returns the insert.
5. Every known restriction enzyme is checked in left and right roles.
6. The curated 19-enzyme picker remains small; the 109-enzyme catalogue is
   for finding sites, not for overwhelming the picker.
7. Isoschizomers may share a result only when site and cut geometry agree;
   neoschizomers must not be merged by recognition sequence alone.
8. Circular vectors, features, reads, and primer read ranges wrap through
   coordinate zero.
9. A vector's cloning cassette may be on the minus strand. pET-21a's NdeI
   and XhoI coordinates cannot be interpreted by ascending number alone.
10. A vector-derived His tag counts only if the tag sequence came from the
    vector portion, not from the insert.
11. Tm is reaction-specific. SSD annealing, PCR annealing portions, and full
    tailed primers are not interchangeable conditions.
12. Reverse-read chromatogram qualities must be mapped back through the
    orientation flip before attaching confidence to differences.

Golden SSD outputs live in `gsynth_engine/tests/test_ssd_golden.py`. PCR and
cloning-primer behaviour lives in `gsynth_engine/tests/test_pcr.py`. Read the
test docstrings as the scientific specification before changing either.

---

## 7. Current implemented scope

| Area | Implemented |
| --- | --- |
| Optimisation | Host codon usage, CAI, GC windows, repeats/homopolymers, restriction-site avoidance, protein-preserving repair |
| SSD | Coding/non-coding paths, NdeI handling, optional stop removal, tag/linker/protease options, duplex and warning metadata |
| Merzoug assembly | PCR-free paired fragments, automatic 4–8 nt overhang widening, uniqueness/palindrome/one-mismatch checks, re-ligation verification |
| PCR | Conventional and cloning modes, adaptive annealing footprints, Tm/Ta, clamps, internal-site blocking, digest simulation, ORF warnings |
| Clone | Vector digest, insert ligation, observed seams, recut check, ORF/tag outcomes, maps and export |
| Check | Ligation calculations, sequencing primers, text-read placement, `.ab1` parsing, Mott trimming, Q-aware differences and peaks |
| Compare | Gotoh affine-gap alignment; global, local and semi-global; DNA/protein; reverse-complement search |
| Import/export | SnapGene `.dna`, GenBank, FASTA; GenBank/FASTA/oligo FASTA/CSV/protocol outputs |
| Accounts/projects | JWT login/refresh/revocation, per-user projects, Django admin |
| Interface | Home, Design, Optimise, PCR, Clone, Check, Compare, Learn, Projects, Viewer, Help |

Only pET-21a(+) and pET-21(+) ship with authoritative vector sequences.
pET-28a(+), pET-22b(+), pET-32a(+), pGEX-4T-1, and pUC19 are catalogued but
must be imported from the lab's or supplier's authoritative file. Never type
or invent a vector sequence to remove that requirement.

Type IIS enzymes are deliberately unsupported by the current cut model.
Their cuts fall outside the recognition sequence and require a model change,
not another row in the existing table.

---

## 8. Workspace persistence and design state

The redesign at `462a20d` is the current production UI. Results and form
state now remain visible while the signed-in user navigates between sections.
Each working page has an explicit **Clear** action. The state provider lives
inside the authenticated shell, so signing out destroys it and prevents data
from crossing accounts.

Persistence is currently **memory-only**:

- survives Design → Optimise/PCR → Design and equivalent client-side route
  navigation;
- does not survive a browser refresh, tab close, browser restart, or sign-out;
- uploaded `File` objects are also memory-only.

If “until I clear it” is intended to include refresh/restart, implement a
versioned `sessionStorage` or `localStorage` layer for serialisable values and
make the privacy choice explicit. Do not attempt to serialise uploaded trace
files or place sensitive sequence data in persistent browser storage without
the user's approval.

Final design QA evidence is on the remote branch:

```text
agent/finalize-design-assets
ddc6fdcdad47afd792b3ba43328df1ae07b46677
```

It adds `design-qa.md` and canonical JPG reference/implementation/comparison
images. The application code itself is already on `main`; that branch is
evidence, not a required production code dependency.

---

## 9. Deployment

`render.yaml` creates:

- `gsynth-api`: Django/Gunicorn web service in Frankfurt.
- `gsynth-app`: static Vite site with React Router fallback.

Production variables:

| Variable | Owner | Purpose |
| --- | --- | --- |
| `DATABASE_URL` | Render secret | Supabase Session Pooler PostgreSQL URL |
| `DJANGO_SECRET_KEY` | Render-generated secret | Django signing key |
| `ALLOWED_HOSTS` | API | Exact Render API hostname |
| `CORS_ALLOWED_ORIGINS` | API | Exact frontend origin |
| `VITE_API_BASE` | frontend build | Absolute API base URL |
| `PYTHON_VERSION` | Render | Python 3.12 |

Never commit actual values. Production intentionally refuses to start when
critical settings are absent or unsafe rather than falling back to ephemeral
SQLite.

The Learn/Ollama feature is optional. The Render blueprint does not deploy an
Ollama service and does not set an external `OLLAMA_BASE_URL`; therefore Learn
should be considered unavailable in production unless an Ollama-compatible
service is intentionally provisioned. The rest of the product is independent
of it. Local defaults are `http://localhost:11434` and model `llama3.1`.

Safe deployment sequence:

1. Make a small branch; never work directly on remote `main`.
2. Run every verification command in section 2.
3. Keep database migrations backwards-compatible with the currently running
   API where possible.
4. Merge to `main`; Render auto-deploys both services.
5. Watch the API build, migrations, boot, and `/api/health/`.
6. Open the frontend, sign in, and exercise the changed workflow.
7. If boot/health fails, roll back to the previous Render deploy and diagnose
   without exposing environment values in logs or chat.

Full click-by-click instructions are in `django_app/DEPLOY.md`.

Any credentials previously pasted into a chat or screenshot must be treated
as compromised and rotated at the provider. This repository and handover do
not contain them.

---

## 10. Repository state at handover

Remote references after an explicit refresh:

```text
origin/main                         462a20de6b3cc0136f5a37c34b41a37b267fff95
origin/agent/finalize-design-assets ddc6fdcdad47afd792b3ba43328df1ae07b46677
```

The local audit workspace was on `agent/finalize-design-assets` at local
commit `bdaacc0a6ed0d2638a4c5fda9ee56af0164d85b0`. Its tracked tree contains
the same design-evidence content as the remote branch, but the commit ID is
different. Use the remote commit above on the workstation.

Five local PNGs were untracked in the audit workspace:

```text
frontend/design-qa/comparison-1.png
frontend/design-qa/comparison-final.png
frontend/design-qa/implementation-1.png
frontend/design-qa/implementation-final.png
frontend/design-reference/design-1.png
```

They are duplicate/intermediate local artifacts. Canonical committed evidence
uses JPG. Do not add the PNGs unless there is a deliberate asset decision.

This revised `HANDOVER.md` is the intended transfer artifact. It is not part
of the production commit named above until explicitly committed and pushed.

---

## 11. Prioritised next work

### P0 — external release evidence

1. Deploy the hardened tree so Render applies `projects.0003`, then repeat the
   authenticated production smoke workflow.
2. Run `docs/USABILITY_STUDY.md` with 5–8 independent bench scientists.
3. Run the physical construct matrix in `docs/WET_LAB_VALIDATION.md` and retain
   gels, raw `.ab1` files, provenance manifests, and reviewer sign-off.

### P1 — remaining compatibility decision

- Decide whether Factor Xa needs exact legacy DNA or only the IEGR peptide;
  record the decision in a golden test before changing the implementation.

### P1 — dependency maintenance

- Keep `package-lock.json` authoritative; use Node 22.22.2 or newer within the
  declared engine range and `npm ci` in CI/deployment.
- Repeat npm and Python advisory audits on every controlled upgrade.

### P1 — production Learn decision

Either provision a supported external Ollama-compatible endpoint with an
explicit privacy/cost assessment, or label Learn as local-only/disabled in
production. Do not send unpublished sequences to a third-party model without
the user's informed approval.

### P2 — optional product work

- Persist serialisable workspace state across refresh only if the user wants
  that stronger meaning of “until Clear”.
- Import authoritative files for the five unbundled vectors.
- Add publication-quality SVG/PDF export for maps and duplexes.
- Consider assembly-list filters for very large fragment sets.

Do not start an interactive sequence editor or Type IIS support without a
clear user requirement; both are substantial scope changes.

---

## 12. Working conventions and known traps

- Use British spelling in product prose and comments.
- Comments explain why; tests state molecular claims.
- User errors say what is wrong and what to do next.
- Severity follows bench consequence: a strategy-destroying internal site is
  a problem; a suboptimal but usable primer is a warning.
- Use one Ruff configuration at repository root.
- Generated enzyme tables and logo assets are changed through their generator,
  not hand-edited.
- Do not group enzymes only by recognition sequence; cut geometry matters.
- Do not quote a terminal end from metadata when it can be measured from the
  molecule.
- In affine-gap traceback, preserve the active matrix layer.
- Use a wide, seeded codon distribution in stress fixtures; short repeating
  generators create artificial k-mer failures.
- Do not apply `zip(..., strict=True)` mechanically; several pairs are
  deliberately offset by one element.
- Visual work requires screenshots at desktop and phone widths plus browser
  console inspection. Passing component tests cannot prove a duplex or map is
  visually correct.
- Preserve unrelated uncommitted work. Stage explicit paths only.

Performance is also a correctness boundary because authenticated design
endpoints accept large sequences. Preserve these shapes:

- codon repair compares candidate changes inside a local window, not against
  the whole gene for every codon;
- banded read alignment iterates and allocates only the band, not a full
  sequence-length row for every base;
- Merzoug overhang selection checks a precomputed exclusion set rather than
  comparing every new junction with every previous junction.

Historical reference timings were about 0.4 s for 3 kb optimisation, 0.01 s
for a 1 kb read, 0.02 s for a 2.4 kb assembly, and 1.3 s for a 200 kb
assembly. Re-benchmark before and after changing those algorithms; do not
treat the old numbers as a current performance certification.

Useful reading order for a new Codex:

1. `HANDOVER.md`
2. `CLAUDE.md`
3. `gsynth_engine/ssd.py`
4. `gsynth_engine/pcr.py`
5. `gsynth_engine/merzoug.py`
6. `gsynth_engine/tests/test_ssd_golden.py`
7. `gsynth_engine/tests/test_pcr.py`
8. `gsynth_engine/tests/test_merzoug.py`
9. `django_app/apps/design/`
10. `frontend/src/state/WorkspaceStateContext.tsx`

When the user says that a scientific result or drawing looks wrong, inspect
the molecule first and the UI second. A plausible rendering can still encode
the wrong sequence.

---

## 13. Definition of done for the next Codex

A change is ready only when:

- its biological choice is explicit;
- engine tests assert the intended molecule and include a case that would
  fail if the protection were removed;
- API output is checked against the engine rather than against a hand-written
  “plausible” value;
- TypeScript, frontend tests, and production build pass;
- affected pages are exercised in a real browser at desktop and phone widths;
- navigation away and back preserves results until Clear where applicable;
- no credential or private sequence was introduced into source, fixtures,
  screenshots, logs, or the commit message;
- `git diff` contains only intended files;
- deployment health is verified after merge.
