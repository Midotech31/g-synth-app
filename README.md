# G-Synth

Software for traceable nucleic-acid design, synthesis ordering, cloning and
post-sequencing validation.

It is not a general sequence editor. Every part of it exists because a step of
that workflow was being done by hand, and the checks it performs are the ones
that, when skipped, cost a fortnight.

## The workflow

```
  gene  →  Optimise  →  Design  →  Hybridise  →  Clone  →  Check
           for the      oligos +   verify the    ligate    ligation ·
           host         protocol   duplex/ends   vector    reads

                         Compare: alignment or physical hybridisation
```

**Optimise** rewrites coding DNA or back-translates a peptide for its expression
host. Automatic start handling preserves a peptide without N-terminal Met as a
mature product, recognises an existing initiator Met, and lets the user override
the biological role when necessary. A complete ORF receives exactly one ATG;
the encoded protein is verified after every rewrite. Fifteen versioned profiles cover six bacterial systems
(*E. coli*, *B. subtilis*, *P. putida*, *L. lactis*, *C. glutamicum* and
*S. coelicolor*), four yeasts (*S. cerevisiae*, *K. phaffii*, *K. lactis* and
*Y. lipolytica*), mammalian human/CHO proxies, insect Sf9/Sf21 and S2 proxies,
and *N. benthamiana*. They are calculated from a committed September 2021
FDA HIVE-CUTs/CoCoPUTs snapshot, using RefSeq when available and a disclosed
GenBank fallback otherwise. The UI shows taxon, dataset, CDS count, codon count,
GC and scope for the selected profile. An experiment-specific highly expressed
reference-gene set can override these species-wide profiles for a strain, cell
line or tissue. The reported profile-relative score is not presented as an
expression-yield prediction. The enzymes you will clone with are also inputs,
because a gene that carries an internal NdeI site cannot be cloned NdeI/XhoI
however favourable its host codons are.

**Design** implements two original workflows. **Small Sequence Design (SSD)**
emits one forward/reverse synthesis pair for a compact construct.
**Extended Sequence Design (ESD)** tiles a longer construct into F/R pairs
with complementary 4–8 nt overhangs that ligate in one order, with no PCR.
Both hand back the exact oligos to order, a bench protocol, and
the hybridisation view: both strands drawn aligned, with the overhangs
showing. Nothing can be downloaded until re-ligating the fragments in silico
reproduces the construct base for base, on both strands — and until the two
outer ends, read off that molecule, are the sticky ends the chosen enzymes
leave, on the right strands.

To our knowledge, G-Synth is the first disclosed automation of SSD and ESD as
these terms are defined here: synthesis-order workflows whose exported
molecules must reconstruct both intended strands exactly before release. This
specific claim does not imply priority over general gene-design, DNA-assembly
or end-to-end construction software.

**Hybridise** is the mandatory molecular gate between SSD/ESD and cloning.
Both order molecules remain entered 5′→3′; G-Synth draws them antiparallel,
distinguishes paired bases, mismatches and exposed 5′/3′ ends, and always shows
both a compact end/core summary and a nucleotide-level double-strand view. A
Design handoff runs this check automatically and transfers only an exact duplex
and its enzyme assignments to cloning. Alignment remains available as a separate
mode for similarity comparison and may introduce gaps; hybridisation never
hides a bulge or cohesive end behind a gap.

**Restriction cloning** cuts a vector and puts the construct in. Its enzyme
assignments can be changed without altering the transferred bases; digestion
must then prove the new end geometry compatible before ligation. pET-21a(+) and pET-21(+)
ship with their sequences; any other backbone is imported — SnapGene `.dna`,
GenBank or FASTA — and checked against the catalogue entry, so pasting
pET-28a while pET-21a is selected is caught rather than cloned into. Each seam
is drawn as the two ends that made it, so "the overhangs match" can be checked
instead of believed. A synchronized Vector / Insert / Product workbench keeps
map, sequence, annotations, enzymes, oligos and the diagnostic gel on one
coordinate system. Product inspection and export remain unavailable until the
reviewed compatible ends are explicitly ligated. The restriction inventory is
searchable, defaults to unique cutters plus the cloning pair, and can expose every
multi-cutter on demand. The coordinate-level Annotated view shows regulatory
features, cassette parts, strand direction and codon-aligned translation,
including expression loci that cross the circular origin. Every feature can
be created, named, edited or deleted; exact matches to a curated library of
common promoters, operators, tags, linkers and cleavage motifs are proposed
for review rather than silently asserted. Saved edits survive into GenBank.

The cloning selectors cover 109 non-redundant cut geometries representing 289
commercial enzyme names. Isoschizomers remain searchable aliases, so the map
does not duplicate identical sites; HindIII is included with its verified
5′-AGCT cohesive end.

**Primer design and PCR simulation** are supporting tools. They distinguish a primer's hybridising 3′ region
from its deliberately unpaired 5′ cloning tail, show how that tail enters the
product after extension, and calculate complete diagnostic-digest fragments.
Predicted agarose gels offer explicit generic 100 bp, 1 kb and broad-range
marker sizes. Every gel is permanently labelled as an in-silico size
prediction, not an experimental image.

**Check** closes the loop: ligation amounts (molar, because at equal mass a
5.4 kb vector outnumbers a 150 bp insert thirty-six to one), sequencing
primers that sit back from the insert rather than at it, and a comparison of
the reads that come back against the design — in either orientation, with a
substitution reported as the residue it changes. Uploaded ABI/AB1 traces are
quality-trimmed and displayed in a reference-aligned chromatogram viewer with
consensus, coordinates, strand direction, base calls, Phred-quality shading,
four-channel peaks and explicit mismatches. Forward and reverse reads are
oriented and assembled before consensus coverage is calculated; raw consensus
coverage, bidirectional overlap, overlap agreement and quality-gated coverage
are reported separately. Only evidence admitted by the selected quality gate
is drawn as confidence-covered.

**Compare** contains both pairwise alignment and physical hybridisation. The
two modes share inputs but keep their distinct scientific interpretations.

## Running it

Everything below was run from a fresh clone before being written down.

### With Docker — one command

```bash
git clone https://github.com/Midotech31/g-synth-app.git
cd g-synth-app/django_app
docker compose up --build
```

Open <http://localhost:5173>. This starts Postgres, Redis, the API and the
frontend. **Learn** is a bundled, searchable molecular-biology knowledge
library, so it works offline without a model download, external service or API
key.

### Without Docker

```bash
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

The frontend proxies `/api` to `:8000`; point it elsewhere with
`VITE_API_TARGET`. Learn contains fixed, reviewable explanations of PCR,
cloning, annotation, expression and post-sequencing validation. It does not
send questions or sequences to an AI service. Optimise, Design, PCR, Clone,
Check and Compare/Hybridise are deterministic calculation tools; Learn is their scoped
reference companion.

### The engine on its own

The engine has no dependencies outside the standard library, so it can be
installed and used from a script or a notebook without any of the above:

```bash
pip install -e .
```

```python
from gsynth_engine import design_extended_sequence

plan = design_extended_sequence(my_gene, enzyme_pair="NdeI / XhoI", is_coding=True)
assert plan.verify() == []          # empty means the oligos re-ligate to the design
```

### Tests

```bash
python -m pytest gsynth_engine/tests -q     # 1,113 — the biology
cd django_app && python -m pytest -q        # 265 — the HTTP layer
cd frontend && npm test                     # 91 — the interface
```

All three run in CI on every push. The engine's suite is the definition of
correctness: it is where the golden examples live, where every one of the 109
enzymes is checked in both positions, and where the property the whole method
rests on is asserted — that the designed fragments re-ligate into the
construct exactly.

Scientific scope and validation materials:

- [`docs/BIOLOGICAL_ASSUMPTIONS.md`](docs/BIOLOGICAL_ASSUMPTIONS.md) — what is checked, and what still requires bench evidence.
- [`docs/WORKED_CLONING_EXAMPLES.md`](docs/WORKED_CLONING_EXAMPLES.md) — cohesive, mixed-polarity, blunt, and blocking cases.
- [`docs/SCIENTIFIC_REFERENCES.md`](docs/SCIENTIFIC_REFERENCES.md) — sources behind the assumptions.
- [`publication_evidence/codon_host_profile_validation.json`](publication_evidence/codon_host_profile_validation.json) — machine-readable host-table provenance, completeness, numerical distinctness and protein-invariance checks.
- [`publication_evidence/peptide_and_enzyme_validation.json`](publication_evidence/peptide_and_enzyme_validation.json) — current glargine A/B protein-to-SSD-to-clone evidence, peptide-start decisions, host-wise translation invariance, enzyme-name coverage and HindIII geometry.
- [`docs/USABILITY_STUDY.md`](docs/USABILITY_STUDY.md) — the human-validation protocol and release criteria.
- [`docs/WET_LAB_VALIDATION.md`](docs/WET_LAB_VALIDATION.md) — physical construct-to-sequencing evidence and acceptance criteria.
- [`docs/ACCESSIBILITY.md`](docs/ACCESSIBILITY.md) — current accessibility and responsive-design guarantees and validation limits.
- [`docs/AI_USE_DISCLOSURE.md`](docs/AI_USE_DISCLOSURE.md) — final, auditable declaration of generative-AI assistance and human responsibility.

## What is in here

| Path | What it is |
|---|---|
| `gsynth_engine/` | The biology. Dependency-free Python, no Django, no HTTP. Every design decision and every check lives here, with its tests. |
| `django_app/` | A thin HTTP layer over the engine, plus accounts and per-user projects. It validates requests and serialises results; it computes nothing. |
| `frontend/` | React + TypeScript workspace. Draws what the engine returns. |

The separation is deliberate: the engine can be imported, tested and trusted
without a web server, and a bug in the biology has one place to be.

## Design notes

A few decisions that are easy to reverse by accident:

- **Restriction enzymes are stored as cut positions**, not as the bases each
  oligo should carry. A pair-specific oligo remainder is correct only for the
  enzyme pair it was derived from and silently produces mismatched duplexes
  when reused for a different enzyme.
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

## Contributing

`CONTRIBUTING.md` describes where a change belongs, what a test here is
expected to assert, and the properties that must not regress.

## Citing this work

If G-Synth contributed to published results, please cite it. `CITATION.cff`
carries the machine-readable record, and GitHub renders it as a **Cite this
repository** button. Cite the version you actually used — the software is
tagged per release, and a construct designed under one version cannot be
checked against another.

## Licence

MIT — see `LICENSE`. Use it, modify it, build on it, including commercially;
the only condition is that the copyright notice travels with it.

## Author

**Prof. Merzoug Mohamed** — Full Professor
Genomics Technology Platform, Higher School of Biological Sciences of Oran
<mohamed.merzoug.essbo@gmail.com>
