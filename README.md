# G-Synth

Software that automates one laboratory's gene synthesis and cloning workflow:
turning a gene into oligos you can order, a plasmid you can build, and the
checks that say it is the construct you designed.

It is not a general sequence editor. Every part of it exists because a step of
that workflow was being done by hand, and the checks it performs are the ones
that, when skipped, cost a fortnight.

## The workflow

```
  gene  →  Optimise  →  Design  →  Clone  →  Check      Compare
           for the      oligos +   into a    ligation ·  two
           host         protocol   vector    primers ·   sequences
                                             reads
```

**Optimise** rewrites a gene for the organism that will express it. The protein
never changes. The enzymes you will clone with are an input, because a gene
that translates beautifully and carries an internal NdeI site cannot be cloned
NdeI/XhoI.

**Design** builds the cassette — sticky ends, ATG, 6×His, linkers, protease
site — and cuts it into oligo pairs for **Merzoug assembly**: fragments
hybridised in F/R pairs with complementary 4–8 nt overhangs, ligated in order,
no PCR at any step. It hands back the oligos to order, a bench protocol, and
the hybridisation view: both strands drawn aligned, with the overhangs
showing. Nothing can be downloaded until re-ligating the fragments in silico
reproduces the construct base for base, on both strands — and until the two
outer ends, read off that molecule, are the sticky ends the chosen enzymes
leave, on the right strands.

**Clone** cuts a vector and puts the construct in. pET-21a(+) and pET-21(+)
ship with their sequences; any other backbone is imported — SnapGene `.dna`,
GenBank or FASTA — and checked against the catalogue entry, so pasting
pET-28a while pET-21a is selected is caught rather than cloned into. Each seam
is drawn as the two ends that made it, so "the overhangs match" can be checked
instead of believed.

**Check** closes the loop: ligation amounts (molar, because at equal mass a
5.4 kb vector outnumbers a 150 bp insert thirty-six to one), sequencing
primers that sit back from the insert rather than at it, and a comparison of
the reads that come back against the design — in either orientation, with a
substitution reported as the residue it changes.

**Compare** aligns two sequences that are not assumed to be the same thing.

## Running it

Everything below was run from a fresh clone before being written down.

```bash
git clone https://github.com/Midotech31/g-synth-app.git
cd g-synth-app

# Backend
python -m venv .venv && source .venv/bin/activate
pip install -r django_app/requirements.txt
cd django_app
python manage.py migrate --settings=config.settings.dev
python manage.py runserver --settings=config.settings.dev      # :8000

# Frontend, in a second terminal
cd frontend
npm install
npm run dev                                                     # :5173
```

Open <http://localhost:5173>, create an account, and the workspace is there.
The frontend proxies `/api` to `:8000`; point it elsewhere with
`VITE_API_TARGET`.

### Tests

```bash
python -m pytest gsynth_engine/tests -q     # 476 — the biology
cd django_app && python -m pytest -q        # 184 — the HTTP layer
```

The engine's suite is the definition of correctness. It is where the golden
examples live, where every one of the 19 enzymes is checked in both positions,
and where the property the whole method rests on is asserted: that the
designed fragments re-ligate into the construct exactly.

## What is in here

| Path | What it is |
|---|---|
| `gsynth_engine/` | The biology. Dependency-free Python, no Django, no HTTP. Every design decision and every check lives here, with its tests. |
| `django_app/` | A thin HTTP layer over the engine, plus accounts and per-user projects. It validates requests and serialises results; it computes nothing. |
| `frontend/` | React + TypeScript workspace. Draws what the engine returns. |
| `app.py`, `modules/`, `utils/` | **The original Streamlit application**, kept as it was. Superseded by the above, not deleted. |

The separation is deliberate: the engine can be imported, tested and trusted
without a web server, and a bug in the biology has one place to be.

## Design notes

A few decisions that are easy to reverse by accident:

- **Restriction enzymes are stored as cut positions**, not as the bases each
  oligo should carry. Storing the latter — as earlier versions did — is
  correct only for the pair it was checked against, and silently produces
  mismatched duplexes for every other enzyme.
- **Melting temperatures come from the nearest-neighbour model** (SantaLucia
  1998) under the conditions of the annealing reaction the protocol
  prescribes, not from base composition and not at a generic primer dilution.
  The two differ by about 7 °C.
- **Vector sequences ship only when they came from an authoritative file.**
  A transcription error in 5 000 bases is invisible and would poison every
  design made against that backbone.
- **Substitution matrices are generated data, not literals.** BLOSUM62 is
  read from Biopython's copy: 576 values typed by hand is 576 chances to be
  quietly wrong.
- **A terminal end is measured, not quoted.** Every terminal value in a plan
  is copied from the design when the fragments are built, so comparing them
  compares a label with itself. The two ends are read back off the assembled
  duplex and checked against the enzyme table, polarity included.
- **Junction overhangs widen when four bases run out.** Each one placed rules
  out every overhang within a base of it, and of its partner, so only 22 can
  coexist — while a 2.4 kb gene needs 26. The design goes to five bases rather
  than failing, and says so; the oligos are the same length either way.

## Deployment

`django_app/DEPLOY.md` covers Supabase + Render, both on free tiers.
