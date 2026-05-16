# G-Synth 3.0

**Advanced gene-synthesis & molecular-biology platform.** Complete ground-up rewrite in layered, testable Python that replaces G-Synth 2025.2.0.

## Highlights

- **`gsynth_core`** — headless library, pip-installable, fully typed, no UI dependencies.
- **Scientifically correct algorithms** where G-Synth 2.x had placeholders:
  - SantaLucia 1998 Tm with Owczarzy salt/Mg/dNTP corrections (validates within 1.3°C of literature M13 primer).
  - Doench 2016 CFD off-target scoring (faithful matrix from Table S19).
  - Multi-objective codon optimization (CAI + GC window + homopolymer + restriction-site avoidance + RBS-mimicry).
  - True 6-frame ORF finder (G-Synth 2.x missed reverse strand).
  - Full IUPAC reverse complement (G-Synth 2.x lost R/Y/S/W/K/M info).
- **80-enzyme curated restriction database** (plus REBASE fallback via Bio.Restriction).
- **Cloning suite** — Gibson (Tm-matched overlaps), Golden Gate (BsaI/BsmBI with CIDAR overhangs), classical restriction.
- **Format support** — FASTA, GenBank (multi-record), SnapGene `.dna`, NCBI Entrez, UniProt (with the corrected REST endpoint).
- **Streamlit UI** demonstrating the architecture.

## Install

```bash
pip install -e .             # core library
pip install -e ".[ui]"       # + Streamlit UI
pip install -e ".[all]"      # everything
```

## Quick start

```python
from gsynth_core.sequence import reverse_complement, find_orfs
from gsynth_core.thermo import melting_temperature
from gsynth_core.codon import optimize, OptimizationParams
from gsynth_core.restriction import find_sites
from gsynth_core.crispr import design_guides, CasType

rc = reverse_complement("ATGCRYSWATG")              # "CATWSRYGCAT" — IUPAC-aware
orfs = find_orfs("ATG...", min_length_aa=50)        # all 6 frames
tm  = melting_temperature("GTAAAACGACGGCCAGT")      # 57.3°C (M13 ~56°C)
opt = optimize("MKRLAVF...",
               params=OptimizationParams(organism="e_coli",
                                         avoid_sites=("GAATTC", "GGATCC")))
guides = design_guides("ATGCTAG...", cas=CasType.SpCas9, min_on_target=0.5)
```

## UI

```bash
streamlit run app.py
```

## Run tests

```bash
pytest
```

Or, without pytest installed:

```bash
python run_tests.py
```

## Architecture

```
gsynth/
├── gsynth_core/           Pure library — no UI deps
│   ├── constants.py       Single source of truth (genetic code, IUPAC, codon tables)
│   ├── sequence/          ops + 6-frame ORF + translation
│   ├── thermo/            SantaLucia 1998 NN Tm
│   ├── codon/             real CAI, multi-objective optimization
│   ├── restriction/       80-enzyme DB + digestion simulator
│   ├── primers/           primer3-py + cloning primers
│   ├── crispr/            Doench 2016 + CFD off-target
│   ├── cloning/           Gibson, Golden Gate, restriction
│   ├── alignment/         NW global, SW local, MSA
│   ├── io/                FASTA, GenBank, SnapGene, NCBI, UniProt
│   └── annotation/        auto-detect plasmid features
├── app.py                 Streamlit entry point (at repo root for Streamlit Cloud)
├── pages/                 Streamlit multi-page navigator
├── gsynth_ui/             UI helpers (state, theme)
├── tests/                 pytest
├── assets/                Logo and brand assets
├── Dockerfile
└── .github/workflows/ci.yml
```

## License

MIT © 2025 Dr. Mohamed Merzoug — Higher School of Biological Sciences of Oran
