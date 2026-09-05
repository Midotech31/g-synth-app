# Wet-lab design-to-sequencing validation

**Status:** protocol ready; no physical run has been claimed.

The automated release gate exercises Design → fragment assembly → Hybridisation → restriction cloning → sequencing-primer design
→ exact read placement against the released plasmid. It proves software
coherence, not that digestion, ligation, transformation, expression, or Sanger
sequencing succeeded at the bench. A physical run requires laboratory staff,
reagents, instruments, samples, and recorded primary evidence.

## Minimum construct matrix

Run at least three independent constructs:

| Construct | Strategy | Purpose |
| --- | --- | --- |
| Lab NdeI/XhoI coding insert | pET-21a(+), NdeI/XhoI | Primary validated expression path |
| Lab mixed-polarity insert | Validated 5′/3′ cohesive pair | Duplex polarity and orientation |
| Internal-site negative control | One extra selected site | Release block and failure interpretation |

Use the laboratory's intended sequences and vector lots. Do not substitute the
synthetic regression fixture for a physical construct.

## Required run record

For each construct retain:

- G-Synth provenance manifest and preflight codes;
- exported oligo/primer order file and cloning worksheet;
- vector identity, lot, concentration, enzyme lots, buffer, incubation, and
  cleanup method;
- uncut, single-cut, double-cut, and ligation-control gel images with ladder;
- colony count and colony-PCR or diagnostic-digest evidence;
- raw ABIF (`.ab1`) or SCF files from both directions (or additional primers until the entire
  intended insert and both junctions are covered);
- G-Synth verification report with coverage, gaps, differences, quality, and
  final five-state verdict;
- reviewer name/date and any deviation from the exported plan.

## Acceptance criteria

1. The ordered molecules and primer sequences equal the exported records.
2. Diagnostic digest band sizes agree with the worksheet within the gel method's
   declared tolerance.
3. Both vector–insert junctions and 100% of the intended insert are covered by
   placed reads.
4. The final state is `fully_verified`; partial coverage is never accepted.
5. No confident sequence difference remains unexplained.
6. The observed plasmid sequence hashes to the released design or the deviation
   is reviewed and a new release record is created.

## Stop conditions

Stop and investigate on an unexpected digest band, an unplaced read, a coverage
gap, a confident difference, mixed chromatogram peaks, a provenance mismatch,
or a discrepancy between ordered and exported molecules. Do not relabel a
partial or ambiguous result as verified.
