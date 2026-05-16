"""Pairwise and multiple sequence alignment."""
from gsynth_core.alignment.pairwise import Alignment, align_global, align_local
from gsynth_core.alignment.msa import MSAResult, multiple_alignment
__all__ = ["Alignment", "align_global", "align_local", "MSAResult", "multiple_alignment"]
