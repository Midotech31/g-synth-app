# Changelog

All notable changes to G-Synth are documented here. Releases follow semantic
versioning.

## [1.0.0] - 2026-08-29

### Added

- Complete design-to-clone workflow with scientific preflight gates.
- Codon optimisation that preserves translation and avoids selected internal
  restriction sites.
- PCR-free Merzoug oligonucleotide assembly with in-silico re-ligation.
- Restriction-enzyme cloning, vector validation, recut checks and export.
- Conventional and cloning-primer design with digest simulation.
- Sanger ABI/AB1 and SCF parsing, quality-aware placement, difference windows,
  and a reference-aligned four-channel chromatogram viewer.
- Reproducibility, wet-lab validation, accessibility, usability and explicit
  AI-use disclosure materials.

### Validation

- 1,056 engine, 232 API and 57 frontend automated tests.
- Insulin-glargine A/B retrospective design and Sanger case study.

### Known limitations

- Archived glargine traces provide partial quality-admitted insert coverage and
  do not establish full-clone verification.
- Experimental validation does not establish expression, folding, bioactivity,
  clinical equivalence or manufacturing economics.
