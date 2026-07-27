# gsynth_engine

The part of G-Synth that is genuinely G-Synth's: the design logic for
synthesising and cloning genes by **oligo hybridisation and ligation**.

Pure Python, no framework, no UI, no dependencies. Both the Streamlit app
and the Django API import it, so there is one implementation of the biology
and one set of tests protecting it.

```
peptide / gene
   │
   ├─ SSD              forward + reverse oligos carrying the exact sticky
   │                   ends of a chosen restriction pair, plus ATG, 6×His,
   │                   flexible linkers and a protease cleavage site
   │
   ├─ Merzoug Assembly long constructs split into oligo pairs joined by
   │                   complementary 4–8 nt overhangs, ligated successively.
   │                   No PCR at any step.
   │
   └─ Protocol         an order-ready oligo sheet and a bench protocol
```

## Use it

```python
from gsynth_engine import design_merzoug_assembly
from gsynth_engine.protocol import bench_protocol, order_sheet_csv

plan = design_merzoug_assembly(
    "GGCATCGTGGAACAGTGCTGCACC...",     # your insert
    enzyme_pair="NdeI / XhoI",
    cleavage_site="Thrombin",
    target_oligo_length=90,
    overhang_length=4,
)

assert plan.verify() == []                 # re-ligates to the design
print(order_sheet_csv(plan, construct_name="pGS-EntA"))
print(bench_protocol(plan))
```

For a short insert that needs no fragmentation, use the SSD stage alone:

```python
from gsynth_engine import design_small_sequence

result = design_small_sequence(insert, enzyme_pair="NdeI / XhoI")
result.forward, result.reverse      # the two oligos
result.left_overhang                # "TA"  — NdeI
result.right_overhang               # "TCGA" — XhoI
result.coding_region                # forward from the ATG onwards
result.segments                     # labelled parts, for display
```

## What the tests protect

`pytest gsynth_engine` — 135 tests.

**Golden tests** (`test_ssd_golden.py`) reproduce the two worked examples in
the G-Synth specification base for base, both strands. They are the
definition of correct. If one fails, the engine would send someone to the
bench with the wrong DNA — fix the code, not the expectation.

They also pin NdeI's subtlety: the forward oligo starts `TATG`, the 5'-TA
overhang and the ATG **overlap by one base**, and ligation restores `CATATG`
because the cut vector supplies the `CA`. That is easy to "simplify" by
accident.

**Duplex integrity** (`test_duplex_integrity.py`) checks, for all 19 enzymes
in both the left and right positions, that the two oligos actually anneal
without a mismatch. G-Synth 2.x stored one pair of post-cut remainders per
enzyme, but what each oligo must carry depends on which end the enzyme sits
at — the stored values were right for NdeI-as-left and XhoI-as-right (the
validated pair) and wrong for every other enzyme, producing duplexes that
would not have annealed. Remainders are now derived from the cut positions,
so all 19 are correct.

**Assembly** (`test_merzoug.py`) checks the method's promise directly:
ligating the oligo pairs in order reproduces the design on both strands, for
inserts from 220 to 1500 bp and overhangs of 4, 6 and 8 nt. Junction
overhangs must be unique, non-palindromic, not homopolymers, not all-AT or
all-GC, and never equal to a terminal enzyme overhang — every one of those
rules exists because breaking it scrambles an assembly silently.

## Design notes

**Fragments are cut from the finished construct, not built up to it.** The
full SSD duplex is designed first, then cut in silico: top strand at the
junction, bottom strand *k* nucleotides further along. The stagger is what
creates the overhangs. Two things then come for free — the terminal ends are
exactly the SSD enzyme overhangs, and re-ligating every fragment must
reproduce the construct, which `verify()` asserts before any plan is
returned.

**Junction placement is a design decision, not arithmetic.** Overhangs
determine assembly order, so junctions start at evenly spaced ideal
positions and are nudged outward until the overhang passes every rule. When
no position within the search window works, the engine raises with the
reason rather than returning a plan that would misassemble.

**Errors are written for the person at the bench**, not the developer:
"The insert contains an internal XhoI site — it will be cut during cloning"
rather than a stack trace.
