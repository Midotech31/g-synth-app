# Biological assumptions and release limits

G-Synth performs deterministic *in-silico* design checks. It does not replace
manufacturer instructions, experimental controls, sequence verification, or a
qualified scientist's review of the intended expression system.

## Restriction cloning

- A vector is releasable only when each selected enzyme cuts it exactly once,
  the two enzymes differ, and the observed insert duplex ends anneal to the
  observed vector ends.
- Cut geometry is derived from the configured enzyme table. The table checksum
  is stored in every provenance record so a design can be tied to the exact
  definitions used to make it.
- Internal sites in a PCR product or insert are blocking because digestion
  would produce more than the intended fragment.
- Site regeneration is a review item, not a universal requirement. Losing a
  site can be intentional, but it changes diagnostic-digest options.
- The app predicts sequence compatibility. It does not predict star activity,
  methylation sensitivity, buffer compatibility, partial digestion, enzyme
  lot performance, or ligation efficiency. Confirm these against current
  supplier documentation.

## Coding constructs

- “Supplies ATG” means the *retained bases after cleavage* place an initiating
  ATG at the insert boundary; an ATG merely visible in a recognition sequence
  is not enough.
- Reading-frame checks use the chosen origin and junction sequence. Protein
  expression still depends on promoter, ribosome-binding/Kozak context, host,
  toxicity, mRNA structure, codon usage, and protein stability.
- Tag and protease-site predictions are sequence predictions, not evidence of
  accessibility, cleavage efficiency, solubility, or purification yield.
- Factor Xa recognition is defined at peptide level as IEGR↓ (with Asp also
  tolerated at P3 by the supplier). G-Synth fixes the DNA spelling to
  `ATC-GAA-GGT-CGT` for reproducible oligo output and E. coli-compatible codon
  use. Other synonymous spellings encode the same peptide, but the fixed
  sequence avoids the rare E. coli AGG codon and makes exported oligos
  reproducible. This exact choice is protected by a golden molecular test.

## Oligo assembly and PCR

- In-silico reassembly must reproduce both designed strands exactly.
- Internal assembly overhangs must be unique and reverse-complement orthogonal
  within the design. This reduces misassembly but does not guarantee yield.
- Primer Tm is a model under stated ionic conditions. Annealing temperature,
  extension time, polymerase choice, and additives require empirical tuning.
- A tailed cloning primer hybridises through its 3′ annealing region in the
  first cycle. Its 5′ clamp and restriction site are intentionally unpaired;
  they become double-stranded only after extension and subsequent cycles.
- Virtual gels plot sequence-derived fragment sizes on a logarithmic migration
  axis. Band intensity, smearing, partial digestion, supercoiled topology and
  matrix/buffer effects are not inferred and require an experimental gel.
- Long-oligo synthesis quality, secondary structure, and supplier-specific
  purification remain experimental constraints.

## Sequencing verification

- “Fully verified” requires complete requested-region coverage and no detected
  differences. Agreement over a partial region is explicitly “partial match.”
- A low-quality chromatogram call is evidence to repeat or inspect, not proof
  that the design base is correct.
- A sequence match does not verify plasmid topology, copy number, contamination,
  phenotype, or expression. Retain appropriate positive, negative, vector-only,
  and no-ligase controls.

## Reproducibility

Every generated project carries the engine version, parameter checksum,
sequence checksum, vector checksum when applicable, enzyme-table checksum, and
generation time. Raw molecular inputs are represented by hashes in the
provenance manifest rather than duplicated there.
