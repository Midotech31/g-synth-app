# Sequence explorer and annotation review — 2026-09-07

## Changes

- The detailed viewer now exposes the entire molecule, coordinate navigation,
  region navigation, 30/60/90 bases per row, vertical resizing and an expanded
  view. Long records use bounded rendering; a regression reaches base 200,000.
- Clone scans both strands automatically. Detected candidates can be inspected
  on the map and sequence and explicitly retained. Saved projects are scanned
  automatically on opening; candidate acceptance remains explicit.
- Shine–Dalgarno proposals require a possible downstream initiation codon in
  the transcript direction. The 4–14 nt spacer is a screening heuristic, not
  a universal biological rule. Overlapping SD spellings are deduplicated;
  motifs inside an annotated same-strand CDS are not proposed as initiation
  sites. Reverse-strand and circular-origin cases are tested.
- The T7 terminator motif follows the transcriptional orientation of the
  existing pET reference. Detected status and its evidence survive annotation
  saving and GenBank export/reimport.
- GenBank regulatory annotations use `regulatory` with `regulatory_class`.
  Qualifier wrapping no longer inserts spaces inside controlled terms.
- Renaming a CDS preserves its translation bounds. Annotation dialogs retain
  keyboard focus, unknown feature types remain editable, unstranded features
  are not labelled forward, and circular selections respect both endpoints.
- Viewer loading and clipboard/download errors no longer leave misleading
  success or stale-project state.

## Scientific references

1. [Estrada et al., 2024, Unraveling the plasticity of translation initiation in prokaryotes](https://pubmed.ncbi.nlm.nih.gov/38206950/).
   The SD interaction occurs between an upstream mRNA region and the anti-SD
   region of 16S rRNA. Initiation architectures vary among organisms; motif
   identity does not establish function or expression efficiency.
2. [INSDC feature table, version 11.4, April 2026](https://www.insdc.org/submitting-standards/feature-table/).
   Defines regulatory features and their qualifiers, strand/location semantics,
   and the distinction between experimentally supported and inferred features.
3. [INSDC regulatory-class vocabulary](https://www.insdc.org/submitting-standards/controlled-vocabulary-regulatoryclass/).
   Provides `ribosome_binding_site`, `promoter` and `terminator` terminology.

On the plasmid DNA map, an SD annotation identifies the DNA region corresponding
 to the motif in the bacterial mRNA. “Upstream” follows the 5′→3′ transcript,
not increasing plasmid coordinates on both strands. An SD motif is not the
entire RBS, and its absence does not rule out bacterial translation.

## Verification scope

The review combines source inspection of the annotation, rendering, import,
export and interaction paths with the complete existing engine, API and
frontend regression suites. Browser regressions exercise design, hybridization,
cloning, primer revalidation, mobile navigation, accessibility and initial
layouts of the scientific workspaces. The CI artifact `browser-evidence`
contains the captured views.

These checks do not establish experimental regulatory activity or exhaustive
correctness of every possible input, dependency, function and visual state.
The motif collection remains deliberately limited to curated common elements;
it is not a genome annotation pipeline or a universal promoter predictor.
