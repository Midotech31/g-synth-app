# utils/alignment_utils.py

"""
Shared alignment utilities for pairwise and multiple sequence alignment.
Wrapped around Biopython pairwise2 and command-line wrappers for MUSCLE/Clustal.
"""

import subprocess
import shutil
import logging
from typing import Tuple, List, Optional

logger = logging.getLogger("G-Synth:ALIGNUTILS")

try:
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    USING_BIOPYTHON_ALIGN = True
except ImportError:
    USING_BIOPYTHON_ALIGN = False
    logger.warning("Biopython pairwise2 not available.")

# Check if `muscle` and `clustalo` are in PATH
HAS_MUSCLE = shutil.which("muscle") is not None
HAS_CLUSTALO = shutil.which("clustalo") is not None

def pairwise_align(
    seq1: str,
    seq2: str,
    method: str = "global",
    match: int = 2,
    mismatch: int = -1,
    gap_open: int = -0.5,
    gap_extend: int = -0.1
) -> List:
    """
    Perform pairwise alignment between two sequences.
    method: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman).
    Returns list of alignment objects.
    """
    if not USING_BIOPYTHON_ALIGN:
        raise ImportError("Biopython pairwise2 not installed.")
    if method == "global":
        return pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend)
    else:
        return pairwise2.align.localms(seq1, seq2, match, mismatch, gap_open, gap_extend)

def run_muscle(input_fasta: str, output_fasta: str) -> bool:
    """
    Run MUSCLE on input_fasta, write result to output_fasta.
    Returns True if successful.
    """
    if not HAS_MUSCLE:
        logger.error("MUSCLE not found in PATH.")
        return False
    try:
        subprocess.run(["muscle", "-in", input_fasta, "-out", output_fasta], check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"MUSCLE error: {e}")
        return False

def run_clustalo(input_fasta: str, output_fasta: str) -> bool:
    """
    Run Clustal Omega on input_fasta, write result to output_fasta.
    """
    if not HAS_CLUSTALO:
        logger.error("Clustal Omega not found in PATH.")
        return False
    try:
        subprocess.run(["clustalo", "-i", input_fasta, "-o", output_fasta, "--force"], check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Clustal Omega error: {e}")
        return False

####################################################################
# End of alignment_utils.py
####################################################################
