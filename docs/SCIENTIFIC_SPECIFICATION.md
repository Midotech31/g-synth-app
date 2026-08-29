# G-Synth scientific specification

**Status:** canonical software specification  
**Scope:** calculations performed by G-Synth; physical success still requires
the wet-lab evidence in `WET_LAB_VALIDATION.md`.

## Authority and change control

When sources disagree, use this order:

1. an approved laboratory SOP or signed bench example;
2. golden molecular tests in `gsynth_engine/tests/`;
3. this specification;
4. interface copy and handover notes;
5. the 2025 executable and the initial prompt.

The initial prompt is historical design input. It contains a correct PCR-free
Merzoug description followed by an incompatible Gibson/PCR description, and a
sentence saying NdeI removes two bases even though all six worked examples
remove the complete template `ATG`. Neither contradiction may override tested
molecular behaviour.

Coordinates are 0-based and half-open in the engine/API. The interface converts
them to 1-based inclusive ranges for bench readability. DNA is written 5′→3′
unless a strand display labels otherwise.

## Small Sequence Design

For non-coding input, the standard cassette is:

```text
left retained enzyme sequence
ATG unless the left enzyme supplies it
GSS linker · 6×His · SSG linker
optional protease site
input sequence
right retained enzyme sequence
```

For coding input, the sequence must begin with `ATG`. If the left retained
enzyme sequence already supplies the start codon (notably NdeI), the complete
template `ATG` is removed so the expressed protein starts with one methionine.
Optional stop removal truncates at the first in-frame TAA, TAG or TGA.

NdeI/XhoI golden oligos are exact compatibility contracts. Other enzyme ends
are derived from cut geometry in their left or right molecular role; the
legacy executable's generic reusable remainder pair is not authoritative.

## Protease sites

Protease choices are explicit DNA sequences. `Factor Xa` is the modern IEGR
coding sequence. `Factor Xa (legacy DNA)` reproduces the 2025 executable's
synonymous IEGR DNA. They are protein-equivalent but never silently
interchanged in an order record.

## Merzoug assembly

Merzoug assembly is PCR-free. Each purchased fragment is a complementary
forward/reverse oligo pair. Internal joins use unique 4–8 nt cohesive ends and
ligase; restriction-derived ends occur only at the outer construct boundaries.
Junctions must not be palindromic, reused, reverse-complement-reused, or within
one mismatch of another accepted junction. Both reassembled strands must equal
the SSD construct before supplier files can be exported.

It is not Gibson assembly, PCR stitching, a fixed 20 bp overlap method, or an
`AGCT`-padding method.

## PCR and cloning primers

Automatic mode chooses 18–30 nt annealing footprints toward 60 °C. Manual mode
accepts exact 15–60 nt footprints and retains all specificity, GC, homopolymer,
hairpin and dimer warnings. Tailed cloning primers contain:

```text
validated terminal clamp + complete recognition site + annealing footprint
```

With NdeI and a template-leading `ATG`, the default annealing footprint begins
at codon two because `CATATG` supplies the initiator. Keeping both ATGs is an
explicit Met–Met opt-in. Tm/Ta are based on the annealing portion under PCR
conditions; whole-oligo Tm is reported separately. The complete product is
simulated, cut, and checked for additional recognition sites. A blocked pair
may only be replaced by an alternative that passes the same full simulation;
vector uniqueness remains a separate Clone check.

## Translation, ORFs and reverse complement

Translation reports +1, +2, +3, −1, −2 and −3 frames. Complete ORFs begin with
ATG and end with TAA, TAG or TGA. Nested in-frame starts remain distinct
candidates. Reverse-strand ORF coordinates are mapped back to the original
top-strand coordinate system. ORF FASTA/CSV exports include strand, frame,
range, DNA and protein. Reverse complement rejects ambiguous bases rather than
guessing them.

## Optimisation, cloning and verification invariants

- Codon optimisation never changes the translated protein.
- Selected restriction sites and configured motifs are repaired synonymously
  or reported unresolved.
- A vector is circular and may carry its cloning cassette on either strand.
- Ligation requires compatible observed ends, correct polarity and orientation.
- Re-cutting the recombinant plasmid with the selected pair returns the insert.
- A design is not sequencing-verified unless the requested region is fully
  covered and contains no unexplained confident difference.
- Tm is reaction-specific; SSD annealing, PCR and sequencing-primer conditions
  are not interchangeable.

## Unsupported or externally validated scope

Type IIS enzymes require cuts outside their recognition sites and are not
represented by the current cut model. Physical digestion, ligation,
transformation, expression and sequencing success are never inferred from an
in-silico pass. Usability claims require independent bench scientists; wet-lab
claims require retained primary evidence and reviewer sign-off.
