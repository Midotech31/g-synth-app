"""
gsynth_engine — the G-Synth design engine.

Pure Python. No web framework and no UI. This package holds the part of
G-Synth that is genuinely G-Synth's: the oligo design logic for synthesising
and cloning genes by hybridisation and ligation.

    peptide / gene
        → SSD          Small Sequence Design: one forward/reverse order pair
                       with the exact sticky ends for a chosen restriction
                       pair, plus optional tag, linkers and cleavage site
        → ESD          Extended Sequence Design: longer constructs split into
                       orderable oligo pairs joined by complementary 4–8 nt
                       overhangs, with no PCR at any step
        → bench        an order-ready oligo list

Everything here is covered by tests, including golden tests that reproduce
the worked examples from the G-Synth specification base for base. Those
tests exist so that no future refactor — in any framework — can silently
change the sequences the lab orders.
"""

__version__ = "1.0.0"

from gsynth_engine.esd import (
    ESDResult,
    OligoPair,
    design_extended_sequence,
)
from gsynth_engine.hybridization import HybridizationResult, Overhang, hybridize
from gsynth_engine.pcr import PcrPrimer, PcrResult, design_pcr
from gsynth_engine.sequence import (
    clean_dna,
    gc_content,
    is_palindrome,
    reverse_complement,
)
from gsynth_engine.ssd import SSDResult, design_small_sequence
from gsynth_engine.verify import ConsensusReport, assemble_consensus

__all__ = [
    "ESDResult",
    "ConsensusReport",
    "HybridizationResult",
    "OligoPair",
    "Overhang",
    "PcrPrimer",
    "PcrResult",
    "SSDResult",
    "__version__",
    "clean_dna",
    "assemble_consensus",
    "design_extended_sequence",
    "design_pcr",
    "design_small_sequence",
    "gc_content",
    "hybridize",
    "is_palindrome",
    "reverse_complement",
]
