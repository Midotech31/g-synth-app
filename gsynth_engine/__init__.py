"""
gsynth_engine — the G-Synth design engine.

Pure Python. No Streamlit, no Django, no UI. This package holds the part of
G-Synth that is genuinely G-Synth's: the oligo design logic for synthesising
and cloning genes by hybridisation and ligation.

    peptide / gene
        → SSD          forward + reverse oligos with the exact sticky ends
                       for a chosen restriction pair, plus His-tag, linkers
                       and a protease cleavage site
        → Merzoug      long constructs split into oligo pairs joined by
          Assembly     complementary 4–8 nt overhangs, ligated successively,
                       with no PCR at any step
        → bench        an order-ready oligo list

Everything here is covered by tests, including golden tests that reproduce
the worked examples from the G-Synth specification base for base. Those
tests exist so that no future refactor — in any framework — can silently
change the sequences the lab orders.
"""

__version__ = "1.0.0"

from gsynth_engine.merzoug import (
    AssemblyPlan,
    OligoPair,
    design_merzoug_assembly,
)
from gsynth_engine.sequence import (
    clean_dna,
    gc_content,
    is_palindrome,
    reverse_complement,
)
from gsynth_engine.ssd import SSDResult, design_small_sequence

__all__ = [
    "AssemblyPlan",
    "OligoPair",
    "SSDResult",
    "__version__",
    "clean_dna",
    "design_merzoug_assembly",
    "design_small_sequence",
    "gc_content",
    "is_palindrome",
    "reverse_complement",
]
