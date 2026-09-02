# Construct Workbench validation

G-Synth uses one synchronized workspace for the vector, insert and ligated product.
Changing molecule or inspection mode never changes the underlying design record:
map, sequence, annotations, restriction sites, synthesis oligos, diagnostic gel and
the selection inspector are projections of the same coordinates.

## State gates

The product is not treated as a molecule until the user has reviewed both cohesive
junctions and selected **Ligate compatible ends**. Before that action, Product,
Gel, Save, GenBank, FASTA and the bench worksheet are unavailable. After ligation,
the product becomes the active circular molecule and the gated evidence becomes
available. This prevents a planned construct from being confused with a simulated
ligation product.

## Restriction-site inventory

Recognition sites are searched on both strands and across the coordinate origin of
circular molecules. Palindromic matches are de-duplicated; non-palindromic motifs
are found through the recognition sequence or its reverse complement. The default
table contains unique cutters and the active cloning pair. Search, cloning-pair-only
and multi-cutter controls provide the complete inventory without presenting
hundreds of low-value rows at once.

For the bundled 5,443 bp pET-21a(+) sequence, the browser and Python engine both
report 732 occurrences across 72 enzyme specifications. HindIII is a unique site at
0-based coordinate 172 (1-based position 173) with the exact `AAGCTT` recognition
sequence.

## Editable product evidence

Product features can be created, named, typed, stranded, recoloured, edited,
removed and restored. Coordinates are validated against the circular product;
origin-crossing spans may wrap once, but cannot exceed the molecule. Reviewed
annotations are included in project saves and GenBank export. CDS translation
coordinates and truncation status survive the API round trip.

## Verification record

- 1,113 biology-engine tests passed.
- 265 Django/API tests passed.
- 92 frontend tests passed.
- 1,470 automated tests passed in total.
- Ruff, TypeScript, migration drift and production frontend build passed.
- npm and Python dependency audits reported no known vulnerabilities.
- Browser exercise covered pre-ligation locks, explicit ligation, cross-view
  selection, HindIII, enzyme filtering, annotation edit/remove/undo, gel display,
  desktop reflow and mobile reflow with no document-level horizontal overflow.
- Browser console reported no warnings or errors during the exercised workflow.
