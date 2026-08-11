# Contributing to G-Synth

This is laboratory software. A wrong number here becomes a wrong oligo
ordered, a plasmid that will not close, or a protein that does not bind the
column — and the person finds out three weeks later. Everything below exists
because that already happened at least once.

## Where a change belongs

```
gsynth_engine/   the biology. Dependency-free Python. No Django, no HTTP.
django_app/      a thin layer that validates, calls the engine, serialises.
frontend/        React. Draws what the engine returns.
```

**No biology outside `gsynth_engine/`.** If you are computing a sequence in a
view, a serializer or a component, it belongs in the engine — that is where it
can be tested without a web server, and where a bug has one place to be.

## Adding an endpoint

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
cd frontend && npm test                     # the interface
```

All three run in CI on every push, and the engine carries a coverage floor.

Write the test as a claim about the molecule, not about the code. Every
docstring in `gsynth_engine/tests/` says what the check is protecting
against; keep that up, because six months later the *why* is the only part
that still matters.

A test that cannot fail is not a test. `test_gates_can_fail.py` exists to
prove the safety checks catch real errors rather than merely passing — if you
add a guard, add the case that trips it.

### Properties that must never regress

- Re-ligating a Merzoug design in silico reproduces the construct base for
  base, **on both strands**. `AssemblyPlan.verify()` runs before any plan is
  returned; nothing downloads until it is empty.
- The assembled fragments present **the sticky ends the chosen enzymes leave**,
  sequence and polarity, at both outer ends — measured off the molecule, not
  read back from the labels the design wrote.
- Codon optimisation never changes the protein.
- Recutting a recombinant plasmid with the same pair returns the insert.
- Every enzyme is checked in **both** positions, left and right.
- The curated 19 are re-checked against REBASE on every run; a hand-typed
  cut position is otherwise invisible when wrong.

## Decisions that are easy to reverse by accident

Each of these was a real defect. Reintroducing one is silent — the arithmetic
still adds up, and the result is still confidently wrong.

**Enzymes are stored as cut positions**, never as the bases an oligo should
carry. What each oligo carries is *derived* per role, because it differs
depending on which end the enzyme sits at.

**Polarity depends on which end, not only which strand.** A protruding top
strand is a 5′ overhang at a fragment's left end and a 3′ overhang at its
right. Deriving it from the strand alone reports NdeI as a 3′ cutter.

**A vector's cassette may read on the minus strand.** In pET-21a NdeI is at
236 and XhoI at 157 — the left-hand enzyme cuts *after* the right-hand one.

**A vector tag counts only where the vector contributed sequence.** Searching
the whole protein for `HHHHHH` finds the insert's own tag and reports it as
the vector's.

**Tm comes from stacking, under the reaction the protocol prescribes** — not
from base composition, and not at a generic primer dilution. Those differ by
about 7 °C.

**Sequences and matrices ship only as verified data.** Vector sequences come
from an authoritative file; BLOSUM62 is generated from Biopython's copy. A
transcription error in 5 000 bases, or in 576 matrix values, is invisible.

**Circular means circular.** Sites, reads, features and primer read-ranges all
wrap past position 0. In every pET construct the insert sits exactly there.

## Performance

Two endpoints were once hangable by any signed-in user. When touching
`codon.py` or `verify.py`:

- Codon repair compares candidates on a **window**, not the whole gene.
- Banded alignment must **iterate only the band** and size its buffers by the
  band.
- In `merzoug.py`, a candidate overhang is checked against a **precomputed
  exclusion set**, not against every overhang already placed.

Measure before and after. 3 kb optimisation ≈ 0.4 s; a 1 kb read ≈ 0.01 s;
a 2.4 kb assembly ≈ 0.02 s.

## Verifying UI work

Screenshots, not assertions. Start both servers, drive the page with
Playwright, and look at the result — three real bugs were found that way and
by no other means. A drawing that is subtly wrong looks like a passing test.

## Style

British spelling in prose. Comments explain *why*, never *what*. Match the
surrounding style rather than importing a new one. `ruff check .` from the
repository root must pass.

Errors reaching a user are written for them: what is wrong, and what to do
about it. `SequenceError` messages are shown verbatim in the interface.

Severity follows consequence. A leftover restriction site blocks the whole
strategy and is a **problem**; a rare codon left in to satisfy a GC window
costs a little translation speed and is a **warning**. Calling both problems
teaches the user to ignore the word.

## Reporting a problem with a design

If G-Synth produced a construct that did not behave at the bench, that is the
most valuable bug report there is. Please include the input sequence, the
enzyme pair, the options set, and what happened — the design is reproducible
from those, and the case becomes a regression test.
