"""gsynth_core — headless bioinformatics library for G-Synth 3.0."""

__version__ = "3.0.1"

from gsynth_core.sequence import (
    reverse_complement, complement, clean_dna, gc_content,
    translate, find_orfs, validate_dna,
)
from gsynth_core.thermo import melting_temperature

__all__ = [
    "__version__",
    "reverse_complement", "complement", "clean_dna", "gc_content",
    "translate", "find_orfs", "validate_dna",
    "melting_temperature",
]
