# Working on G-Synth

This is laboratory software. A wrong number here becomes a wrong oligo
ordered, a plasmid that will not close, or a protein that does not bind the
column — and the person finds out three weeks later. Everything below exists
because that already happened at least once.

## The shape

```
gsynth_engine/   the biology. Dependency-free Python. No Django, no HTTP.
django_app/      a thin layer that validates, calls the engine, serialises.
frontend/        React. Draws what the engine returns.
app.py, modules/, utils/   the original Streamlit app. Superseded, kept.
```

**No biology outside `gsynth_engine/`.** If you find yourself computing a
sequence in a view, a serializer or a component, it belongs in the engine —
that is where it can be tested without a web server, and where a bug has one
place to be.

Adding an endpoint:

1. Logic in the engine, with tests. Callable and verifiable without Django.
2. A serializer with **bounds on every field**. A field with no maximum is a
   denial of service waiting for someone to paste an operon.
3. A view that calls the engine and turns `SequenceError` into a 400. Give it
   `throttle_scope = "design"` if it runs an algorithm rather than a lookup.
4. An HTTP test asserting the response **matches what the engine returned** —
   not that it looks plausible.

## Tests are the specification

```bash
python -m pytest gsynth_engine/tests -q     # the biology
cd django_app && python -m pytest -q        # the HTTP layer
```

Write the test as a claim about the molecule, not about the code. Every
docstring in `gsynth_engine/tests/` says what the check is protecting
against; keep that up, because six months later the *why* is the only part
that still matters.

The load-bearing properties, which must never regress:

- Re-ligating a Merzoug design in silico reproduces the construct base for
  base, **on both strands**. `AssemblyPlan.verify()` runs before any plan is
  returned; nothing downloads until it is empty.
- Codon optimisation never changes the protein.
- Recutting a recombinant plasmid with the same pair returns the insert.
- Every enzyme is checked in **both** positions, left and right.

## Decisions that are easy to reverse by accident

Each of these was a real defect. Reintroducing one is silent.

**Enzymes are stored as cut positions**, never as the bases an oligo should
carry. What each oligo carries is *derived* per role, because it differs
depending on which end the enzyme sits at. Storing one pair per enzyme is
correct only for the pair it was checked against.

**Polarity depends on which end, not only which strand.** The two strands run
in opposite directions: a protruding top strand is a 5' overhang at a
fragment's left end and a 3' overhang at its right. Deriving it from the
strand alone reports NdeI as a 3' cutter.

**A vector's cassette may read on the minus strand.** In pET-21a NdeI is at
236 and XhoI at 157 — the left-hand enzyme cuts *after* the right-hand one.
Assuming otherwise keeps the cloning stuffer and discards the origin and the
marker, with the arithmetic still adding up.

**A vector tag counts only where the vector contributed sequence.** Searching
the whole protein for `HHHHHH` finds the insert's own tag and reports it as
the vector's C-terminal one: a construct that will not bind the column,
described as one that will.

**Tm comes from stacking, under the reaction the protocol prescribes.** Not
from base composition, and not at a generic primer dilution — those differ by
about 7 °C, and composition cannot tell two oligos of equal GC apart.

**Sequences and matrices ship only as verified data.** Vector sequences come
from an authoritative file (a supplier's or the lab's own SnapGene/GenBank
export) and are validated against their catalogue entry. BLOSUM62 is
generated from Biopython's copy. A transcription error in 5 000 bases, or in
576 matrix values, is invisible.

**Affine gaps need the traceback to remember which layer it is *in*,** not
which layer won at each cell. The latter fragments one twelve-base deletion
into four.

**Four-base overhangs run out.** Placing one junction rules out its whole
one-mismatch neighbourhood *and* its partner's, so exactly 22 mutually-usable
4 nt overhangs exist — whichever ones you pick. A 2.4 kb gene at 90 nt oligos
needs 26, and the design used to fail outright at junction 19. It now widens
the overhang, which costs nothing (fragments are cut from a fixed construct,
so the oligos stay the same length) and says so in a warning, because the
number on the form is no longer the number in the tubes. `OVERHANG_SUPPLY` is
measured, not estimated; a test re-derives it from the rules themselves.

**Circular means circular.** Sites, reads, features and primer read-ranges
all wrap past position 0. Clamping at the ends silently truncates them, and
in every pET construct the insert sits exactly there.

## Performance

Two endpoints were once hangable by any signed-in user, by accident as
easily as on purpose. When touching `codon.py` or `verify.py`:

- Codon repair compares candidates on a **window**, not the whole gene — all
  candidates differ in the same three bases, so the rest cancels.
- Banded alignment must **iterate only the band** and size its buffers by the
  band. Allocating a full row per base turns a linear algorithm quadratic.

The same applies in `merzoug.py`: a candidate overhang is checked against a
**precomputed exclusion set**, not against every overhang already placed.
Comparing pairwise is quadratic in junctions — invisible on a peptide, half a
minute on the 200 kb the endpoint accepts.

Measure before and after. 3 kb optimisation ≈ 0.4 s; a 1 kb read ≈ 0.01 s;
a 2.4 kb assembly ≈ 0.02 s.

## Verifying UI work

Screenshots, not assertions. Start both servers, drive the page with
Playwright (`/opt/pw-browsers/chromium`), and look at the result — three real
bugs were found that way and by no other means. A drawing that is subtly
wrong looks like a passing test.

## Conventions

British spelling in prose. Comments explain *why*, never *what*. Match the
surrounding style rather than importing a new one.

Errors reaching a user are written for them: what is wrong, and what to do
about it. `SequenceError` messages are shown verbatim in the interface.

Severity follows consequence. A leftover restriction site blocks the whole
strategy and is a **problem**; a rare codon left in to satisfy a GC window
costs a little translation speed and is a **warning**. Calling both problems
teaches the user to ignore the word.
