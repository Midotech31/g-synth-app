# G-Synth — handover

Everything an agent needs to continue this work well: what the application
is *for*, how it thinks, what is finished, what is not, and the mistakes this
repository has already sprung on somebody.

Read this first. `CLAUDE.md` is the short operational version of the same
rules and is loaded automatically; this is the reasoning behind it.

**Status at handover:** 856 engine tests, 193 HTTP tests, ruff clean,
frontend builds and typechecks. Branch `claude/update-and-fix-bugs-gQoQw`,
head `5378f81`, nothing unpushed.

---

# PART I — WHAT IT IS FOR

## 1. The objective

G-Synth automates **one laboratory's** gene synthesis and cloning workflow.
It is not a general sequence editor and must not become one. Every feature
exists because a step was being done by hand; every check exists because
skipping it once cost somebody a fortnight.

The lab: **Mohamed Merzoug's group**, working on bacteriocins and enterocins
from *Enterococcus*, expressed in *E. coli* from pET vectors. That decides
the defaults — pET-21a(+), NdeI/XhoI, C-terminal 6×His, thrombin — and it
decides what is worth building. A feature that does not serve that workflow
is scope creep even when it is a good idea.

The stated standard, in the user's own words:

> *"Ma logique / la qualité de rendu et accuracy comparable à SnapGene et
> Geneious."*

Three separate obligations, and they are not equally met:

| Obligation | Meaning | State |
|---|---|---|
| **His logic** | Merzoug assembly, his cassette, his vectors — end to end | **Met.** Nothing else does this at all. |
| **Accuracy** | Right by the standards a commercial tool is held to | **Met**, and exceeded in two places (below). |
| **Rendering** | The user *sees* a thing is right, not is told | **~85%.** Maps, duplexes, junctions, traces exist; interactive editing does not. |

**Where G-Synth beats the commercial tools**, and why that must be
preserved:

1. **It refuses to hand over a design it has not verified.**
   `AssemblyPlan.verify()` re-ligates the design in silico and blocks the
   download unless the fragments reproduce the construct base for base on
   both strands. SnapGene will happily export a construct it never
   simulated.
2. **It judges a chromatogram, it does not merely draw one.** SnapGene shows
   you the peak; G-Synth says *"Q7 — check the trace"* and separates a
   mutation from a bad call. Identical letters, opposite conclusions.

---

# PART II — THE SCIENCE, IN ENOUGH DETAIL TO CHANGE CODE SAFELY

## 2.1 Merzoug Assembly

The lab's own method for building a gene from oligos. **No PCR at any step.**

> Each fragment is a short duplex made by annealing one forward and one
> reverse oligo. Adjacent fragments carry complementary 4–8 nt overhangs, so
> they ligate in exactly one order. Fragments are joined pairwise with ligase
> until the full-length insert exists, which is then cloned with the terminal
> restriction pair.

### How the fragments are derived — the key idea

The full-length duplex is built **first** (the SSD cassette), then cut *in
silico*: the top strand at a junction, the bottom strand k nucleotides
further along.

```
top     5'──────────────┐ ┌──────────────3'
bottom  3'──────────┐   └─┘          ┌───5'
                    └──── k nt ──────┘
                    the overhang each junction presents
```

Deriving fragments by **cutting a fixed construct** — rather than building
fragments and hoping they compose — buys two guarantees for free:

1. the terminal ends are *exactly* the SSD enzyme overhangs, and
2. re-ligating every fragment reproduces the construct base for base.

**Never invert this.** If you find yourself constructing oligos first and
assembling second, both guarantees are gone and no test will tell you.

### Why junction placement is not arbitrary

Overhangs decide the assembly order. Two junctions sharing an overhang — or
an overhang that is its own reverse complement — let fragments ligate in the
wrong order or to themselves, and the error is invisible until sequencing.

Junctions are therefore nudged from their ideal position until every
overhang is unique, non-palindromic, distinct from the vector's own sticky
ends, and **not within one mismatch of any other** (NEB's ligase-fidelity
data: overhangs differing at a single position still join measurably).

### The overhang supply is finite and small

Placing one overhang rules out its own one-mismatch neighbourhood *and* its
reverse complement's. At 4 nt that exhausts the alphabet after **22
junctions**, whichever ones you pick. A 2.4 kb gene at 90 nt oligos needs 26.

The design widens the overhang automatically — 5 nt gives 92, 6 nt gives 482
— and reports it as a **warning**, because the number on the form is no
longer the number in the tubes. `OVERHANG_SUPPLY` in `merzoug.py` is
measured, not estimated; a test re-derives it from the rules themselves.

Widening costs nothing at the bench: fragments are cut from a fixed
construct, so a wider stagger moves *where* the cuts fall, not how much has
to be synthesised. Keep that true if you touch the fragment builder.

## 2.2 SSD — Small Sequence Design

The cassette the insert is wrapped in, before fragmentation:

```
[left sticky end][ATG][6×His][linker][protease site][INSERT][linker][right sticky end]
```

Every element is optional and order matters. `ssd.py` is the heart of the
application; `test_ssd_golden.py` pins known-good outputs. Changing a
constant in `constants.py` changes the oligos the lab orders.

## 2.3 Restriction geometry — the part that is easy to get wrong

Enzymes are stored as **cut positions** (`recognition`, `cut_top`,
`cut_bottom`), never as the bases an oligo should carry. What each oligo
carries is *derived per role*, because it differs depending on which end the
enzyme sits at — `left_remainders()` and `right_remainders()`. Storing one
pair of strings per enzyme is correct only for the pair it was checked
against and silently wrong for every other.

**Polarity depends on the side, not only the strand.** The two strands run
in opposite directions, so a protruding top strand is a **5' overhang at a
fragment's left end** and a **3' overhang at its right**. Deriving polarity
from the strand alone reports NdeI as a 3' cutter.

**Two enzyme sets, for two different questions:**

- `RESTRICTION_ENZYMES` — **19**, curated. What a user *picks a cloning pair
  from*. This lab's freezer. Never widen it to the size of the other; a
  hundred-name dropdown is worse than nineteen.
- `ALL_ENZYMES` — **109**, generated from REBASE. What answers *"what else
  cuts here"*. pET-21a has 60 single-cutters; a map drawn from nineteen
  showed 14, and a diagnostic digest could be planned against sites nobody
  was shown.

**Isoschizomers collapse; neoschizomers do not.** NheI and BmtI both cut
`GCTAGC` — NheI leaves `G^CTAG_C` (5'), BmtI leaves `G_CTAG^C` (3').
Grouping by recognition sequence merges them and returns the wrong sticky
end for whichever name loses. Group by *(site, cut_top, cut_bottom)*.

## 2.4 Cloning

`cloning.py` cuts the vector, inserts the construct, closes the plasmid.

**A vector's cassette may read on the minus strand.** In pET-21a, NdeI is at
236 and XhoI at 157 — the left-hand enzyme cuts *after* the right-hand one.
Assuming otherwise keeps the 78 bp cloning stuffer and throws away the
origin and the resistance marker, with all the arithmetic still adding up.
The fix: flip the vector when the smaller arc demands it, mirror the
annotations.

**A vector tag counts only where the vector contributed sequence.**
Searching the whole protein for `HHHHHH` finds the *insert's own* tag and
reports it as the vector's C-terminal one — describing a construct that will
not bind the column as one that will.

## 2.5 Verification against sequencing

`verify.py` places each read against the design (k-mer anchors vote on
orientation and offset, then a **banded** alignment runs only in the located
window) and reports what differs, where, and whether it changes the protein.

`chromatogram.py` reads the `.ab1` trace directly — ABIF is a published
format and the engine takes no runtime dependencies. Having the quality
values changes three things:

1. **Trimming is by quality** (Mott's algorithm), not by a fixed count. A
   fixed 30 discards 60 good bases on a clean read and keeps 30 useless ones
   on a bad one.
2. **Every difference carries the confidence of the base that produced it.**
   Below Q20 it is marked unsupported — at that confidence the basecaller is
   choosing between a peak and its neighbour's shoulder.
3. **The peaks around each difference are returned and drawn**, so the
   judgement is made by looking rather than taken on trust.

A read supplied as letters reports `quality: null`, **not** "fine". An
unmarked difference must mean *unknown*, or the whole distinction is worth
nothing.

---

# PART III — ARCHITECTURE

## 3.1 The rule above all others

```
gsynth_engine/   the biology. Dependency-free Python. No Django, no HTTP.
django_app/      a thin layer that validates, calls the engine, serialises.
frontend/        React + TypeScript. Draws what the engine returns.
tools/           generators for committed assets (logo, enzyme table).
app.py, modules/, utils/   the original Streamlit app. Superseded, kept.
```

**No biology outside `gsynth_engine/`.** If you find yourself computing a
sequence in a view, a serializer or a component, it belongs in the engine —
that is where it can be tested without a web server and where a bug has one
place to be.

The engine has **no runtime dependencies**. Biopython appears only as a
*build-time* source for generated data (BLOSUM62, the enzyme table) and in
the Django layer's file parsing. Do not add an engine dependency; if a
format needs parsing, parse it — `chromatogram.py` is ~200 lines for ABIF.

## 3.2 Adding an endpoint

1. Logic in the engine, with tests. Callable and verifiable without Django.
2. A serializer with **bounds on every field**. A field with no maximum is a
   denial of service waiting for someone to paste an operon.
3. A view that calls the engine and turns `SequenceError` into a 400. Give
   it `throttle_scope = "design"` if it runs an algorithm, not a lookup.
4. An HTTP test asserting the response **matches what the engine returned** —
   not that it looks plausible.

## 3.3 Module map

| Module | Lines | Responsibility | Watch out for |
|---|---|---|---|
| `sequence.py` | 79 | primitives, `SequenceError` | keep dependency-free |
| `constants.py` | 137 | curated 19, `ALL_ENZYMES`, 7 protease sites | cut positions, never strings |
| `enzyme_table.py` | 126 | generated REBASE table (109) | regenerate, never hand-edit |
| `thermo.py` | 224 | SantaLucia 1998 NN Tm, salt corrections | `ANNEALING` ≠ `PRIMER`, ~7 °C apart |
| `ssd.py` | 303 | the cassette | golden tests pin the output |
| `merzoug.py` | 680 | fragmentation, junction placement | `OVERHANG_SUPPLY`, `terminal_ends`, `verify()` |
| `duplex.py` | 398 | hybridisation view, junction views | polarity by side |
| `cloning.py` | 865 | cut, insert, close, validate | minus-strand cassettes; observed vs declared ends |
| `vectors.py` | 493 | 7-vector catalogue, validation | only authoritative sequences ship |
| `codon.py` | 674 | optimisation, CAI (Sharp & Li 1987) | windowed cost — §5 |
| `verify.py` | 500 | read placement, differences | banded DP — §5 |
| `chromatogram.py` | 338 | ABIF, quality, Mott trimming | reverse-read index mapping |
| `align.py` | 425 | Gotoh affine-gap pairwise | traceback layer — §4 |
| `primers.py` | 316 | sequencing primers | read ranges wrap the origin |
| `ligation.py` | 195 | molar ligation arithmetic | molar, not mass |
| `protocol.py` | 261 | order sheet, bench protocol | |
| `genbank.py` | 179 | GenBank / FASTA export | LOCUS line is column-positional |

**HTTP:** 17 design endpoints under `/api/design/`, 7 accounts, 1 sequences.

**Frontend** (~4 560 lines): 7 workspace pages — Design, Optimise, Clone,
Verify, Align, Dashboard, Viewer — plus Login and Signup; 6 components —
DuplexView, JunctionDuplex, TraceView, InsertForm, Logo, Icon.

---

# PART IV — INVARIANTS THAT MUST NEVER REGRESS

Each is enforced by a test. If one breaks, the application is producing
**wrong molecules**, not merely misbehaving.

1. **Re-ligating a Merzoug design in silico reproduces the construct base
   for base, on both strands.** `verify()` runs before any plan is returned;
   nothing downloads until it is empty.
2. **The assembled fragments present the sticky ends the chosen enzymes
   leave** — sequence *and* polarity, at both outer ends, measured off the
   molecule by `terminal_ends`, not read back from the labels the design
   wrote.
3. **Codon optimisation never changes the protein.**
4. **Recutting a recombinant plasmid with the same pair returns the insert.**
5. **Every enzyme is checked in both positions, left and right** — all 109.
6. **The curated 19 are re-checked against REBASE on every run.** A
   hand-typed cut position is otherwise invisible when wrong.
7. **Circular means circular.** Sites, reads, features and primer
   read-ranges wrap past position 0. Clamping silently truncates them, and
   in every pET construct the insert sits exactly there.

## Decisions that are easy to reverse by accident

Each was a **real defect**. Reintroducing one is silent — the code runs, the
tests you thought to write pass, and the molecule is wrong.

- **A terminal end is measured, never quoted.** Every terminal value in a
  plan is copied from the SSD when fragments are built, so a check that
  compares them agrees with itself whatever the oligos spell. Same reason
  `cloning.py` has `_observed_insert_ends`.
- **Affine gaps need the traceback to remember which layer it is *in*,** not
  which layer won at each cell. The latter fragments one twelve-base
  deletion into four.
- **Sequences and matrices ship only as verified data.** Vector sequences
  come from an authoritative file and are validated against their catalogue
  entry. BLOSUM62 is generated from Biopython. A transcription error in
  5 000 bases, or in 576 matrix values, is invisible.
- **Tm comes from stacking, under the reaction the protocol prescribes.**
  Not from base composition, and not at a generic primer dilution —
  composition cannot tell two oligos of equal GC apart.
- **Reverse reads map back through the flip.** The trace runs in the read's
  own direction and the comparison does not. Mapping a difference back
  without undoing the flip reads the quality of the base the same distance
  from the *other* end: a real number, from the wrong base, entirely
  plausible.

---

# PART V — PERFORMANCE ENVELOPE

Two endpoints were once hangable by any signed-in user. The serializers
accept up to 200 kb, so every algorithm must be sane at that size.

| Fix | Before | After |
|---|---|---|
| Codon repair compares candidates on a **window**, not the whole gene | 6.7 s @ 3 kb | 0.42 s |
| Banded alignment **iterates only the band** and sizes buffers by the band | >2 min | 1.9 s |
| Junction placement uses a **precomputed exclusion set**, not pairwise | ~25 s @ 200 kb | 1.3 s |

Reference timings: 3 kb optimisation ≈ 0.4 s · 1 kb read ≈ 0.01 s ·
2.4 kb assembly ≈ 0.02 s · 200 kb assembly ≈ 1.3 s.

**Measure before and after.** A perf test pins the placement bound.

---

# PART VI — WHAT IS DONE

Every item below is implemented, tested, and verified in a browser.

## Workflow

| Area | What works | Where |
|---|---|---|
| **Optimise** | Codon optimisation for the host, CAI (Sharp & Li), GC windows, homopolymer and repeat limits, restriction-site avoidance, rare-codon repair. Protein never changes. | `codon.py`, `/optimise/` |
| **Design** | SSD cassette; Merzoug fragmentation to oligo pairs; automatic overhang widening with a warning; order sheet (CSV); bench protocol; hybridisation view; per-junction views. Gene scale proven to 200 kb. | `ssd.py`, `merzoug.py`, `protocol.py`, `duplex.py` |
| **Clone** | Cut vector, insert, close; minus-strand cassettes handled; 6 validation checks; ORF and tag outcomes; junction seams drawn; restriction map (28 features on pET-21a). | `cloning.py`, `/clone/` |
| **Check** | Sequencing verification in either orientation; `.ab1` chromatogram reading with quality-aware trimming and per-difference confidence; sequencing primers; molar ligation arithmetic. | `verify.py`, `chromatogram.py`, `primers.py`, `ligation.py` |
| **Compare** | Gotoh affine-gap pairwise, global / local / semi-global, DNA and protein (BLOSUM62), reverse-complement search. | `align.py`, `/align/` |
| **Projects** | Save a design or plasmid; reopen with its content; circular / linear / both plasmid maps with features. | `Dashboard.tsx`, `Viewer.tsx` |
| **Export** | GenBank (features preserved), FASTA, oligo FASTA, order-sheet CSV, bench protocol. | `genbank.py`, `protocol.py` |
| **Import** | SnapGene `.dna`, GenBank, FASTA — validated against the catalogue entry, so pasting pET-28a while pET-21a is selected is caught. | `apps/sequences/parsing.py`, `vectors.py` |

## Accounts and safety

JWT with `ver`-claim revocation and blacklist; per-user project isolation;
`ScopedRateThrottle` on 12 algorithmic endpoints; bounds on every serializer
field; 19 security tests.

## Brand and interface

Vector logo generated by `tools/generate_logo.py` in two cuts (detailed, and
a compact one for ≤40 px and the favicon); 9 drawn SVG icons on one 24 px
grid. **No emoji anywhere** — an emoji is a different picture on every
operating system, and sometimes an empty box. Scientific arrows in prose
(`5'→3'`) are typography and stay as text.

## Test inventory

**Engine — 856 tests**

`align` 25 · `chromatogram` 29 · `cloning` 53 · `codon` 35 ·
`duplex_integrity` 5 · `duplex_view` 38 · `enzymes` 9 (parametrised over 109
→ ~570) · `genbank` 24 · `ligation` 17 · `merzoug` 39 · `primers` 17 ·
`protocol` 15 · `ssd_golden` 7 · `thermo` 23 · `vectors` 34 · `verify` 18

**HTTP — 193 tests**

`design` 127 · `sequences` 24 · `security` 19 · `projects` 12 ·
`accounts` 10 · `health` 1

---

# PART VII — WHAT REMAINS

Ranked by what would actually cost the user bench time.

### 1. Five vectors need importing per use — **BLOCKED ON THE USER**

Only pET-21a(+) and pET-21(+) ship with sequences. pET-28a(+), pET-22b(+),
pET-32a(+), pGEX-4T-1 and pUC19 are catalogued — name, length, resistance,
tags, notes — but have no sequence, so the user must import a file each
time.

**Do not invent these sequences.** A transcription error in 5 000 bases is
invisible and poisons every design made against that backbone. They must
come from the user's own SnapGene export or a supplier's file. snapgene.com,
addgene.org and ncbi.nlm.nih.gov are all blocked at the container's egress
policy — do not try to route around it.

*When the files arrive:* drop each in `gsynth_engine/vector_data/<key>.json`
in the same shape as `pET-21a.json`, set `bundled=` on the `VectorSpec`, and
let `validate()` check it against the catalogue entry. One pass.

### 2. Type IIS enzymes — needs a model change, not a table change

BsaI, BsmBI, Esp3I and their kin cut **outside** their recognition sequence.
A cut is currently stored as an offset *into* the site, so `site[cut_top:]`
returns an empty string rather than raising — **wrong ends, not absent
ones**, which is worse. `tools/generate_enzymes.py` excludes them on
purpose.

Supporting them means letting `cut_top`/`cut_bottom` exceed the site length
and teaching the remainder functions to reach into the flanking sequence.
Worth doing only if the lab moves to Golden Gate.

### 3. Interactive sequence editing — deliberate deferral

Everything renders read-only. You cannot click a base, select a region, or
annotate in place. Building a competent editor is months of work and the
user has SnapGene for browsing. **Do not start this without being asked.**

### 4. Smaller, genuinely useful

- **Publication-quality figure export** (SVG/PDF of the map and duplex).
  The renderings exist; only the export does not.
- **Filters on multi-fragment assemblies** — by overhang type, or to only
  the junctions with a warning. Offered to the user, never requested. At 55
  fragments the junction list is long.
- **A custom map surface.** The circular/linear maps come from SeqViz, a
  library. It works; it gives less control than drawing it.

### 5. Deployment — **do not raise unprompted**

`django_app/DEPLOY.md` covers Supabase + Render on free tiers and is
written. The user has said explicitly he is not handing the app to his team
until it is finished. He will raise it.

---

# PART VIII — HOW TO VERIFY YOUR WORK

```bash
python -m ruff check .                      # one config, at the repo root
python -m pytest gsynth_engine/tests -q     # 856 — the biology
cd django_app && python -m pytest -q        # 193 — the HTTP layer
cd frontend && npm run typecheck && npm run build
```

All four must be green before you touch anything, and again before you
commit.

Write each test as a **claim about the molecule**, not about the code. Every
docstring in `gsynth_engine/tests/` says what the check protects against;
keep that up, because six months later the *why* is the only part that still
matters.

**Prove a check can fail.** An earlier compatibility check passed for every
input because it compared a label with itself. When you add an invariant,
add the test that breaks it deliberately.

## UI work: screenshots, not assertions

Start both servers, drive the page with Playwright
(`executablePath: '/opt/pw-browsers/chromium'`), and **look at the result**.
Several real bugs were found that way and by no other means — a drawing that
is subtly wrong looks exactly like a passing test.

```bash
cd django_app && python manage.py runserver --settings=config.settings.dev
cd frontend && npm run dev        # proxies /api to :8000
```

SQLite by default; no database to install. On Windows the venv activation is
`.venv\Scripts\activate`, not `source`.

## Regenerating committed assets

```bash
python tools/generate_logo.py                                  # → mark.svg
python tools/generate_enzymes.py > gsynth_engine/enzyme_table.py
```

Both need Biopython. **Adjust the generator and re-run it**; never hand-edit
the output. That is how vector art and data tables quietly rot.

---

# PART IX — CONVENTIONS

- **British spelling** in prose. Comments explain *why*, never *what*.
- Match the surrounding style rather than importing a new one.
- Errors reaching a user are written **for them**: what is wrong, and what
  to do about it. `SequenceError` messages are shown verbatim in the
  interface.
- **Severity follows consequence.** A leftover restriction site blocks the
  whole strategy and is a **problem**; a rare codon left in to satisfy a GC
  window costs a little translation speed and is a **warning**. Calling both
  problems teaches the user to ignore the word.
- Commit messages say *why*, with the defect named. The history is part of
  the documentation.

---

# PART X — TRAPS THIS REPOSITORY HAS ALREADY SPRUNG

Learned the hard way. Every one of these cost real time.

- **A test filler that repeats.** A codon list cycled by a short-period
  generator gives a sequence with almost no distinct k-mers, and junction
  placement legitimately fails. That is your generator, not a bug. Use
  `random.Random(seed)` over a wide codon set.
- **Blanket `strict=True` on `zip`.** Most zips here are deliberately
  unequal — `zip(fragments, fragments[1:])` is the pairwise idiom. Applying
  it everywhere broke 97 tests. Decide per site.
- **`from __future__ import annotations` hides undefined names.** An
  un-imported type in a signature never raises at runtime and every test
  passes. Only the linter catches it — which is why there is exactly **one**
  ruff config, at the root. A second copy in a subdirectory means
  `ruff check .` from the root silently uses defaults.
- **Diagnostics that scale wrongly.** A "this sequence is repetitive" check
  compared distinct k-mers against the region's *length*; since no sequence
  has more than 4^k distinct k-mers, it labelled every gene past 131 kb
  repetitive. Compare against what the algorithm consumes, not against size.
- **The viewBox is not the artwork.** An untrimmed SVG viewBox is why a logo
  looks small beside its own wordmark and hangs off-centre.
- **Grouping enzymes by recognition sequence.** It merges neoschizomers and
  returns the wrong sticky end. Group by *(site, cut_top, cut_bottom)*.
- **Don't add biology to a view** to save a round trip. It will be wrong in
  a place that has no tests.
- **Don't trust a fixture from one instrument.** `test_chromatogram.py`
  *writes* ABIF files to the published spec rather than shipping one,
  because a fixture from a particular sequencer tests that sequencer.

---

# PART XI — WHERE TO START

1. Read `CLAUDE.md` (the short rules) and this file.
2. Run all four verification commands. Green before you touch anything.
3. Read `gsynth_engine/merzoug.py` and `gsynth_engine/ssd.py` end to end.
   Everything else serves those two.
4. Read `gsynth_engine/tests/test_merzoug.py` — it is the specification of
   the method, written as claims about molecules.

**On working with this user.** He is a working molecular biologist and he
will tell you what is wrong with a design faster than any test will. When he
says a rendering is misaligned or a number looks off, **check the molecule
before checking the code** — twice in this project he was right and the
first instinct was wrong. He writes in French and English; the codebase is
British English and should stay that way.
