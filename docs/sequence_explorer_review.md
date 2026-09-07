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

## SD context and terminator correction — 7 September 2026

An additional source and sequence review found two distinct annotation problems.
An SD-like spelling near an incidental start codon was typed as an RBS without
an associated CDS. It is now an unassigned `misc_feature`; only a motif upstream
of an annotated translation start is proposed as an RBS, still inferred and
functionally unconfirmed. The evidence gives the strand, spacer, start coordinate
and CDS association. The screen remains limited to two spellings and a 4–14 nt
spacer; absence of a match does not exclude an RBS or translation. Sequence-only
vector imports are reassessed against the actual insert start after cloning.

The bundled pET-21 and pET-21a terminator annotations had the opposite direction
to the asymmetric terminator sequence and expression cassette. Their strand
metadata now agree with the existing motif detector: reverse in the source
records, forward in the displayed recombinant. No DNA bases were changed.
Correcting the metadata also eliminates the opposite-strand duplicate proposal.
Previously saved or imported annotations are not silently rewritten.

The 5,490 bp regression retains its RBS at 5,353–5,358 and ATG at 5,367,
with an 8 nt spacer. Its motif near base 142 is unassigned, not the insert RBS.
Short labels now begin with `?`; full candidate status remains in tooltips and
accessible names, and applies in sequence and map views.

Additional references checked:

- [Barendt et al., 2013](https://pubmed.ncbi.nlm.nih.gov/23427812/):
  experimental evidence for context-dependent, non-SD RBS function supports
  distinguishing an SD motif from the broader RBS concept.
- [Calvopina-Chavez et al., 2022](https://academic.oup.com/g3journal/article/12/6/jkac070/6555000):
  the native T7 terminator comprises an RNA stem-loop followed by a 3′ U-rich
  tract. This asymmetry supports the transcriptional orientation correction;
  a terminator annotation does not promise complete termination.
- [INSDC inference qualifiers](https://www.insdc.org/submitting-standards/inference-qualifiers/):
  distinguish non-experimental annotation evidence from functional validation.
