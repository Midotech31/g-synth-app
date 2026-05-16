# Changelog

## 3.0.1 — 2026

Post-release maintenance pass on top of 3.0.0:

- **Restriction database**: filled in BbsI, DraIII, SfiI so the curated set actually contains the documented 80 enzymes (previous count was 77).
- **CRISPR SpCas9-NG scoring**: `_extract_context` no longer hardcodes `pam_len=3`. For SpCas9-NG (2-nt PAM) the 30-mer Doench window now uses the correct PAM length and a compensating trailing tail, instead of silently misaligning the PAM in the scoring matrix.
- **Dockerfile**: source is copied before `pip install` so setuptools' `packages.find` actually discovers `gsynth_core` and `gsynth_ui`. Without this, `pip install ".[ui]"` installed only the dependencies and left the gsynth packages out of site-packages.
- **SnapGene topology**: replaced a substring-matching fallback that could misclassify a linear plasmid as circular when the `dna` dict contained the word "circular" in any field.
- **UI primers page**: renamed local variable `re` (was shadowing the `re` regex module).
- **GenBank writer**: simplified a tautological `value.startswith(value[0])` condition.
- **Streamlit sidebar**: optional logo loading from `assets/logo.png` for branding continuity with G-Synth 2.x.

## 3.0.0 — 2025

Complete ground-up rewrite replacing G-Synth 2025.2.0.

### Scientific corrections (vs G-Synth 2.x)

- **6-frame ORF finder**: G-Synth 2.x's `find_orfs` searched only the 3 forward frames, missing ~50% of ORFs in any plasmid > a few kb. Now searches all 6 frames with correct reverse-strand coordinate conversion.
- **IUPAC reverse complement**: G-Synth 2.x preserved R/Y/S/W/K/M unchanged, producing incorrect complements for any IUPAC-containing sequence. Now uses a full IUPAC complement table.
- **Tm calculation**: G-Synth 2.x averaged Breslauer 1986, SantaLucia 1996, and Sugimoto 1996 NN parameters — three models with incompatible conventions. Now uses SantaLucia 1998 unified parameters with Owczarzy 2008 salt/Mg/dNTP corrections (validated against M13 universal primer literature value).
- **CRISPR scoring**: G-Synth 2.x's `score_spcas9_guide` returned `50 + (len(guide) % 10) * 2` — not biology. Now implements Doench 2016 Rule Set 2 heuristic on-target scoring and faithful CFD (Table S19) off-target scoring.
- **UniProt fetch**: G-Synth 2.x used `?query=sequence:{seq}` which is not a valid UniProt filter and always returned 0 results. Now uses the entry endpoint for accessions and the search API for free-text queries.
- **Codon optimization**: G-Synth 2.x had two "methods" (`cai` and `most_frequent`) that did exactly the same thing — picked the most frequent codon. Now implements a real multi-objective optimizer: CAI + GC window + homopolymer + direct-repeat + restriction-site avoidance + RBS mimicry.
- **Restriction database**: 22 hard-coded enzymes → 80 curated enzymes (including type-IIS for Golden Gate) with Bio.Restriction REBASE fallback (~700 enzymes when Biopython is installed).
- **`shutil` import bug**: G-Synth 2.x's docking module called `shutil.which("vina")` without importing `shutil` — guaranteed NameError on every click.
- **Constants duplication**: `GENETIC_CODE` was defined in three different modules in G-Synth 2.x with risk of divergence. Now centralized in `gsynth_core/constants.py`.

### New features

- **Cloning suite** — Gibson Assembly (Tm-matched overlaps), Golden Gate / MoClo (BsaI/BsmBI with CIDAR orthogonal overhangs, with verification that enzyme sites are removed from the final product), classical restriction + ligation with directional verification.
- **Pairwise alignment** — Gotoh affine-gap implementation of Needleman-Wunsch and Smith-Waterman.
- **MSA** — wrappers for MAFFT / MUSCLE / Clustal Omega when available, with a pure-Python star-alignment fallback.
- **GenBank parser** — pure-Python, multi-record, handles compound locations (`join(...)`, `complement(...)`).
- **SnapGene `.dna` reader** — via optional `snapgene-reader` package.
- **NCBI Entrez fetch** — proper `efetch` integration.
- **Plasmid feature auto-annotation** — 40-feature curated database (promoters, ORIs, resistance markers, affinity tags, protease sites) with fuzzy matching.
- **Streamlit UI** — centralized session state (sequences persist across reruns), unified theme (no per-module CSS).
- **Test suite** — 103 unit tests, all passing.
- **Docker + CI** — production Dockerfile and GitHub Actions workflow.

### Architecture

- Clean separation between `gsynth_core` (headless library) and `gsynth_ui` (Streamlit wrapper).
- Library is fully typed, has no Streamlit dependency, and is pip-installable.
- Optional dependencies (primer3-py, Biopython, snapgene-reader, MAFFT/MUSCLE) — all functionality degrades gracefully when absent.
