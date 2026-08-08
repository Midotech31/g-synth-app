# G-Synth — what it is for, and how it thinks

A handover document. Read this before changing anything; then read
`CLAUDE.md`, which is the short operational version of the same rules.

---

## 1. The objective

G-Synth automates **one laboratory's** gene synthesis and cloning workflow.
It is not a general sequence editor, and it should not become one. Every
feature exists because a step was being done by hand, and every check exists
because skipping it once cost somebody a fortnight.

The lab: Mohamed Merzoug's group, working on **bacteriocins and enterocins
from *Enterococcus***, expressed in *E. coli* from pET vectors. That context
decides the defaults — pET-21a(+), NdeI/XhoI, C-terminal 6×His, thrombin —
and it decides what is worth building. A feature that does not serve that
workflow is scope creep even if it is a good idea.

The stated standard is **"la qualité de rendu et l'accuracy comparable à
SnapGene et Geneious"**. Interpret that as two separate obligations:

- **Accuracy** — the numbers must be right by the same standards a commercial
  tool is held to: real thermodynamics, real cut geometry, real alignment
  algorithms. Not approximations that look plausible.
- **Rendering** — the user must be able to *see* that a thing is right, not be
  told. Where SnapGene draws a plasmid, G-Synth must draw the equivalent.

There is a third obligation the commercial tools do **not** meet, and it is
G-Synth's actual advantage: **it refuses to hand over a design it has not
verified.** SnapGene will happily export a construct it never simulated.
`AssemblyPlan.verify()` runs before any plan is returned, and nothing
downloads until it comes back empty.

---

## 2. The science, in enough detail to change the code safely

### 2.1 Merzoug Assembly

The lab's own method for building a gene from oligos. **No PCR at any step.**

```
Each fragment is a short duplex made by annealing one forward and one
reverse oligo. Adjacent fragments carry complementary 4–8 nt overhangs, so
they ligate in exactly one order. Fragments are joined pairwise with ligase
until the full-length insert exists, which is then cloned with the terminal
restriction pair.
```

**How the fragments are derived — this is the key idea.** The full-length
duplex is built *first* (the SSD cassette), then cut *in silico*: the top
strand at a junction, the bottom strand k nucleotides further along.

```
top     5'──────────────┐ ┌──────────────3'
bottom  3'──────────┐   └─┘          ┌───5'
                    └──── k nt ──────┘
                    the overhang each junction presents
```

Deriving fragments by cutting a fixed construct — rather than by building
fragments and hoping they compose — buys two guarantees for free:

1. the terminal ends are *exactly* the SSD enzyme overhangs, and
2. re-ligating every fragment reproduces the construct base for base.

Never invert this. If you ever find yourself constructing oligos first and
assembling second, you have lost both guarantees.

**Why junction placement is not arbitrary.** Overhangs decide the assembly
order. Two junctions sharing an overhang — or an overhang that is its own
reverse complement — let fragments ligate in the wrong order or to
themselves. The error is invisible until sequencing. So junctions are nudged
away from their ideal position until every overhang is unique,
non-palindromic, distinct from the vector's sticky ends, and not within one
mismatch of any other (NEB's ligase-fidelity data: overhangs differing at a
single position still join at a measurable rate).

**The overhang supply is finite and small.** Placing one overhang rules out
its own one-mismatch neighbourhood *and* its reverse complement's. At 4 nt
that exhausts the alphabet after **22 junctions**, whichever ones you pick.
A 2.4 kb gene at 90 nt oligos needs 26. The design therefore widens the
overhang automatically (5 nt → 92, 6 nt → 482) and says so in a warning.
`OVERHANG_SUPPLY` in `merzoug.py` is measured, not estimated; a test
re-derives it from the rules themselves.

Widening costs nothing at the bench: fragments are cut from a fixed
construct, so a wider stagger moves *where* the cuts fall, not how much has
to be synthesised. Verify this stays true if you touch the fragment builder.

### 2.2 SSD — Small Sequence Design

The cassette the insert is wrapped in, before fragmentation:

```
[left sticky end][ATG][6×His][linker][protease site][INSERT][linker][right sticky end]
```

Every element is optional and order matters. `ssd.py` is the heart of the
application; `test_ssd_golden.py` pins known-good outputs.

### 2.3 Restriction geometry — the part that is easy to get wrong

Enzymes are stored as **cut positions** (`recognition`, `cut_top`,
`cut_bottom`), never as the bases an oligo should carry. What each oligo
carries is *derived per role*, because it differs depending on which end the
enzyme sits at. `left_remainders()` and `right_remainders()` do that
derivation. Storing one pair of strings per enzyme is correct only for the
one pair it was checked against, and silently wrong for every other.

**Polarity depends on the side, not only the strand.** The two strands run in
opposite directions, so a protruding top strand is a **5' overhang at a
fragment's left end** and a **3' overhang at its right**. Deriving polarity
from the strand alone reports NdeI as a 3' cutter. 109 cut specifications are
supported and every one is tested in **both** positions.

### 2.4 Cloning

`cloning.py` cuts the vector, inserts the construct, and closes the plasmid.

**A vector's cassette may read on the minus strand.** In pET-21a, NdeI is at
236 and XhoI at 157 — the left-hand enzyme cuts *after* the right-hand one.
Assuming otherwise keeps the 78 bp cloning stuffer and throws away the origin
and the resistance marker, with all the arithmetic still adding up. The fix:
flip the vector when the smaller arc demands it, and mirror the annotations.

**A vector tag counts only where the vector contributed sequence.** Searching
the whole protein for `HHHHHH` finds the *insert's own* tag and reports it as
the vector's C-terminal one — describing a construct that will not bind the
column as one that will.

### 2.5 Verification against sequencing

`verify.py` places each read against the design (k-mer anchors vote on
orientation and offset, then a **banded** alignment runs in the located
window) and reports what differs, where, and whether it changes the protein.

`chromatogram.py` reads the `.ab1` trace itself — ABIF is a published format
and the engine takes no runtime dependencies. Having the quality values
changes three things:

1. **Trimming is by quality** (Mott's algorithm), not by a fixed count. A
   fixed 30 discards 60 good bases on a clean read and keeps 30 useless ones
   on a bad one.
2. **Every difference carries the confidence of the base that produced it.**
   Below Q20 it is marked unsupported — at that confidence the basecaller is
   choosing between a peak and its neighbour's shoulder. *Identical letters,
   opposite conclusions.* This is the single most valuable thing the module
   does.
3. **The peaks around each difference are returned and drawn**, so the
   judgement is made by looking rather than taken on trust.

A read supplied as letters reports `quality: null`, **not** "fine". An
unmarked difference must mean *unknown*, or the distinction is worthless.

---

## 3. Architecture — one rule above all others

```
gsynth_engine/   the biology. Dependency-free Python. No Django, no HTTP.
django_app/      a thin layer that validates, calls the engine, serialises.
frontend/        React + TypeScript. Draws what the engine returns.
tools/           generators for committed assets (the logo).
app.py, modules/, utils/   the original Streamlit app. Superseded, kept.
```

**No biology outside `gsynth_engine/`.** If you find yourself computing a
sequence in a view, a serializer or a component, it belongs in the engine —
that is where it can be tested without a web server, and where a bug has one
place to be.

The engine has **no runtime dependencies**. Biopython appears only as a
*build-time* source for generated data (BLOSUM62) and in the Django layer's
file parsing. Do not add an engine dependency; if a format needs parsing,
parse it (see `chromatogram.py`, ~200 lines for ABIF).

### Adding an endpoint

1. Logic in the engine, with tests. Callable and verifiable without Django.
2. A serializer with **bounds on every field**. A field with no maximum is a
   denial of service waiting for someone to paste an operon.
3. A view that calls the engine and turns `SequenceError` into a 400. Give it
   `throttle_scope = "design"` if it runs an algorithm rather than a lookup.
4. An HTTP test asserting the response **matches what the engine returned** —
   not that it looks plausible.

---

## 4. Module map

| Module | Lines | Responsibility | Watch out for |
|---|---|---|---|
| `sequence.py` | 79 | primitives, `SequenceError` | keep dependency-free |
| `constants.py` | 133 | curated 19 + `ALL_ENZYMES` (109), 7 protease sites | enzymes as positions, never as strings |
| `enzyme_table.py` | — | generated REBASE table | regenerate, never hand-edit |
| `thermo.py` | 224 | SantaLucia 1998 NN Tm, salt corrections | `ANNEALING` ≠ `PRIMER` — ~7 °C apart |
| `ssd.py` | 299 | the cassette | golden tests pin the output |
| `merzoug.py` | 680 | fragmentation, junction placement | `OVERHANG_SUPPLY`, `terminal_ends`, `verify()` |
| `duplex.py` | 398 | hybridisation view, junction views | polarity by side |
| `cloning.py` | 865 | cut, insert, close, validate | minus-strand cassettes; observed vs declared ends |
| `vectors.py` | 493 | 7-vector catalogue, validation | only authoritative sequences ship |
| `codon.py` | 674 | optimisation, CAI | windowed cost — see §6 |
| `verify.py` | 500 | read placement + differences | banded DP — see §6 |
| `chromatogram.py` | 338 | ABIF reading, quality, trimming | reverse-read index mapping |
| `align.py` | 425 | Gotoh affine-gap pairwise | traceback layer — see §5 |
| `primers.py` | 316 | sequencing primers | read ranges wrap the origin |
| `ligation.py` | 195 | molar ligation arithmetic | molar, not mass |
| `protocol.py` | 261 | order sheet + bench protocol | |
| `genbank.py` | 179 | GenBank/FASTA export | LOCUS line is column-positional |

**17 design endpoints** under `/api/design/`, plus accounts (7) and sequences (1).

**Frontend**: 7 workspace pages (Design, Optimise, Clone, Verify, Align,
Dashboard, Viewer) plus Login and Signup, and 6 components (DuplexView,
JunctionDuplex, TraceView, InsertForm, Logo, Icon).

---

## 5. Load-bearing properties — must never regress

Each is enforced by a test. If one of these breaks, the application is
producing wrong molecules, not merely misbehaving.

1. **Re-ligating a Merzoug design in silico reproduces the construct base for
   base, on both strands.** `AssemblyPlan.verify()` runs before any plan is
   returned; nothing downloads until it is empty.
2. **The assembled fragments present the sticky ends the chosen enzymes
   leave**, sequence *and* polarity, at both outer ends — measured off the
   molecule by `terminal_ends`, not read back from the labels the design
   wrote.
3. **Codon optimisation never changes the protein.**
4. **Recutting a recombinant plasmid with the same pair returns the insert.**
5. **Every enzyme is checked in both positions, left and right.**
6. **Circular means circular.** Sites, reads, features and primer read-ranges
   all wrap past position 0. Clamping silently truncates them, and in every
   pET construct the insert sits exactly there.

### Decisions that are easy to reverse by accident

Each was a real defect. Reintroducing one is **silent** — the code runs, the
tests you thought to write pass, and the molecule is wrong.

- **A terminal end is measured, never quoted.** Every terminal value in a plan
  is copied from the SSD when the fragments are built, so a check that
  compares them agrees with itself whatever the oligos spell. Same reason
  `cloning.py` has `_observed_insert_ends`.
- **Affine gaps need the traceback to remember which layer it is *in*,** not
  which layer won at each cell. The latter fragments one twelve-base deletion
  into four.
- **Two enzyme sets, for two questions.** `RESTRICTION_ENZYMES` is what a
  user *picks a cloning pair from* — nineteen, what this lab keeps.
  `ALL_ENZYMES` is what answers *"what else cuts here"* — 109, because
  pET-21a has 60 single-cutters and a map drawn from nineteen showed 14.
  Never widen the first to the size of the second: a hundred-name dropdown
  is worse than nineteen.
- **Isoschizomers collapse; neoschizomers do not.** NheI and BmtI share
  `GCTAGC` and cut it opposite ways — 5' against 3'. Grouping by recognition
  sequence merges them and hands back the wrong sticky end for whichever
  name loses. Group by *(site, cut_top, cut_bottom)*.
- **Sequences and matrices ship only as verified data.** Vector sequences come
  from an authoritative file (a supplier's or the lab's own SnapGene/GenBank
  export) and are validated against their catalogue entry. BLOSUM62 is
  generated from Biopython's copy. A transcription error in 5 000 bases, or in
  576 matrix values, is invisible.
- **Tm comes from stacking, under the reaction the protocol prescribes.** Not
  from base composition, and not at a generic primer dilution — composition
  cannot tell two oligos of equal GC apart.
- **Reverse reads map back through the flip.** The trace runs in the read's
  own direction and the comparison does not. Mapping a difference back
  without undoing the flip reads the quality of the base the same distance
  from the *other* end: a real number, from the wrong base, entirely
  plausible.

---

## 6. Performance envelope

Two endpoints were once hangable by any signed-in user, by accident as easily
as on purpose. The serializers accept up to 200 kb, so every algorithm must
be sane at that size.

- **Codon repair** compares candidates on a **window**, not the whole gene —
  all candidates differ in the same three bases, so the rest cancels.
  (6.7 s → 0.42 s on 3 kb.)
- **Banded alignment** must iterate *only the band* and size its buffers by
  the band. Allocating a full row per base turns a linear algorithm
  quadratic. (>2 min → 1.9 s.)
- **Junction placement** checks a candidate against a **precomputed exclusion
  set**, not against every overhang already placed. Pairwise is quadratic in
  junctions — invisible on a peptide, half a minute on 200 kb.
  (~25 s → 1.3 s at the endpoint's own maximum.)

Reference timings: 3 kb optimisation ≈ 0.4 s · 1 kb read ≈ 0.01 s ·
2.4 kb assembly ≈ 0.02 s · 200 kb assembly ≈ 1.3 s.

**Measure before and after.** A perf test pins the placement bound.

---

## 7. How to verify your work

```bash
python -m ruff check .                      # one config, at the repo root
python -m pytest gsynth_engine/tests -q     # 856 — the biology
cd django_app && python -m pytest -q        # 193 — the HTTP layer
cd frontend && npm run typecheck && npm run build
```

Write each test as a **claim about the molecule**, not about the code. Every
docstring in `gsynth_engine/tests/` says what the check protects against;
keep that up, because six months later the *why* is the only part that still
matters.

**Prove a check can fail.** An earlier compatibility check passed for every
input because it compared a label with itself. When you add an invariant, add
the test that breaks it deliberately.

### UI work: screenshots, not assertions

Start both servers, drive the page with Playwright
(`executablePath: '/opt/pw-browsers/chromium'`), and **look at the result**.
Several real bugs were found that way and by no other means. A drawing that
is subtly wrong looks exactly like a passing test.

```bash
cd django_app && python manage.py runserver --settings=config.settings.dev
cd frontend && npm run dev        # proxies /api to :8000
```

SQLite by default; no database to install.

---

## 8. Conventions

- **British spelling** in prose. Comments explain *why*, never *what*.
- Match the surrounding style rather than importing a new one.
- Errors reaching a user are written **for them**: what is wrong, and what to
  do about it. `SequenceError` messages are shown verbatim in the interface.
- **Severity follows consequence.** A leftover restriction site blocks the
  whole strategy and is a **problem**; a rare codon left in to satisfy a GC
  window costs a little translation speed and is a **warning**. Calling both
  problems teaches the user to ignore the word.
- **No emoji anywhere in the interface.** An emoji is a different picture on
  every operating system — and sometimes an empty box. All icons are drawn
  SVG on one 24 px grid, one stroke weight, `currentColor` (`Icon.tsx`).
  Scientific arrows in prose (`5'→3'`) are typography and stay as text.
- Generated assets ship **with their generator** (`tools/generate_logo.py`),
  so shapes are adjusted and regenerated rather than hand-edited as path
  data — which is how vector art quietly rots.

---

## 9. State of the application

**Complete and verified:** Optimise · Design (Merzoug + SSD) · Clone ·
Check (sequencing verification, chromatograms, primers, ligation) · Compare ·
Projects · export to GenBank/FASTA/CSV · circular and linear plasmid maps.

**Known gaps, in the order they are worth closing:**

1. **Five of seven vectors need importing per use.** Only pET-21a(+) and
   pET-21(+) ship with sequences. pET-28a, pET-22b, pET-32a, pGEX-4T-1 and
   pUC19 are catalogued but need the user's own SnapGene export each time.
   *Blocked on the user supplying `.dna` files — snapgene.com is not
   reachable from the build container. Do not invent these sequences.*
2. ~~19 enzymes~~ **Done.** 109 distinct cut specifications, generated from
   REBASE by `tools/generate_enzymes.py` and committed. Type IIS enzymes
   (BsaI, BsmBI, Esp3I) are still excluded: they cut *outside* their
   recognition sequence, and the data model stores a cut as an offset into
   the site, so `site[cut_top:]` comes back empty rather than raising.
   Supporting them means extending the model, not extending the table.
3. **No interactive sequence editing.** Everything renders read-only. This is
   a deliberate deferral: building a competent editor is months of work and
   the user has SnapGene for browsing. Do not start it without being asked.

**Deployment** (`django_app/DEPLOY.md`, Supabase + Render) is written but the
user has said explicitly he is not handing the app to his team yet. **Do not
raise deployment unprompted.**

---

## 10. Traps a new agent will fall into

Learned the hard way in this repository.

- **A test filler that repeats.** A codon list cycled by a short-period
  generator produces a sequence with almost no distinct k-mers, and junction
  placement legitimately fails. That is your generator, not a bug. Use
  `random.Random(seed)` over a wide codon set.
- **Blanket `strict=True` on `zip`.** Most zips here are deliberately unequal
  (`zip(fragments, fragments[1:])` is the pairwise idiom). Applying it
  everywhere broke 97 tests. Decide per site.
- **`from __future__ import annotations` hides undefined names.** An
  un-imported type in a signature never raises at runtime and every test
  passes. Only the linter catches it — which is why there is exactly one
  ruff config, at the root. A second copy in a subdirectory means
  `ruff check .` from the root silently uses defaults.
- **Diagnostics that scale wrongly.** A "this sequence is repetitive" check
  compared distinct k-mers against the region's *length*; since no sequence
  has more than 4^k distinct k-mers, it labelled every gene past 131 kb
  repetitive. Compare against what the algorithm consumes, not against size.
- **The viewBox is not the artwork.** An untrimmed SVG viewBox is why a logo
  looks small beside its own wordmark. Measure the bounding box.
- **Don't add biology to a view** to save a round trip. It will be wrong in
  a place with no tests.

---

## 11. Where to start

1. Read `CLAUDE.md` (the short rules) and this file.
2. Run all four verification commands above; they should be green before you
   touch anything.
3. Read `gsynth_engine/merzoug.py` and `gsynth_engine/ssd.py` end to end.
   Everything else serves those two.
4. Read `gsynth_engine/tests/test_merzoug.py` — it is the specification of
   the method, written as claims about molecules.

The user is a working molecular biologist. He will tell you what is wrong
with a design faster than any test will. When he says a rendering is
misaligned or a number looks off, check the molecule before checking the
code.
