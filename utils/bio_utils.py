# utils/bio_utils.py

"""
Shared Bioinformatics Utilities for G-Synth.
This file centralizes core sequence-based functions originally from G-Synth_2025_5_0.py.
"""

import re
import math
import logging
from typing import Tuple, List, Dict, Optional

logger = logging.getLogger("G-Synth:BIOUTILS")

# ────────────────────────────────────────────────────────────────────────────────
# 1. SSD_RESTRICTION_ENZYMES (used by Primer Generator for Cloning)
# ────────────────────────────────────────────────────────────────────────────────
SSD_RESTRICTION_ENZYMES: Dict[str, Dict[str, str]] = {
    "NdeI":    {"recognition": "CATATG",   "cut_forward": "TATG",   "cut_reverse": "CA"},
    "XhoI":    {"recognition": "CTCGAG",   "cut_forward": "C",      "cut_reverse": "TCGAG"},
    "EcoRI":   {"recognition": "GAATTC",   "cut_forward": "AATT",   "cut_reverse": "G"},
    "BamHI":   {"recognition": "GGATCC",   "cut_forward": "GATC",   "cut_reverse": "C"},
    "HindIII": {"recognition": "AAGCTT",   "cut_forward": "AGCT",   "cut_reverse": "T"},
    "SalI":    {"recognition": "GTCGAC",   "cut_forward": "TCGA",   "cut_reverse": "C"},
    "XbaI":    {"recognition": "TCTAGA",   "cut_forward": "CTAG",   "cut_reverse": "A"},
    "NcoI":    {"recognition": "CCATGG",   "cut_forward": "CATG",   "cut_reverse": "G"},
    "KpnI":    {"recognition": "GGTACC",   "cut_forward": "",       "cut_reverse": "GGTAC"},
    "SacI":    {"recognition": "GAGCTC",   "cut_forward": "",       "cut_reverse": "GAGCT"},
    "NotI":    {"recognition": "GCGGCCGC", "cut_forward": "GGCC",   "cut_reverse": "GC"},
    "SpeI":    {"recognition": "ACTAGT",   "cut_forward": "CTAG",   "cut_reverse": "T"},
    "PstI":    {"recognition": "CTGCAG",   "cut_forward": "",       "cut_reverse": "CTGCA"},
    "BglII":   {"recognition": "AGATCT",   "cut_forward": "GATC",   "cut_reverse": "T"},
    "SmaI":    {"recognition": "CCCGGG",   "cut_forward": "",       "cut_reverse": ""},
    "ApaI":    {"recognition": "GGGCCC",   "cut_forward": "",       "cut_reverse": ""},
    "MluI":    {"recognition": "ACGCGT",   "cut_forward": "",       "cut_reverse": ""},
    "EcoRV":   {"recognition": "GATATC",   "cut_forward": "",       "cut_reverse": ""},
    "HpaII":   {"recognition": "CCGG",     "cut_forward": "",       "cut_reverse": ""},
    "SspI":    {"recognition": "AATATT",   "cut_forward": "",       "cut_reverse": ""},
    "DdeI":    {"recognition": "CTNAG",    "cut_forward": "",       "cut_reverse": ""},
    "Bsu36I":  {"recognition": "CCTNAGG",  "cut_forward": "",       "cut_reverse": ""}
}

# ────────────────────────────────────────────────────────────────────────────────
# 2. Clean / Validate DNA Sequence
# ────────────────────────────────────────────────────────────────────────────────
def clean_dna_sequence(seq: str, keep_ambiguous: bool = False) -> str:
    """
    Remove invalid characters from a DNA sequence.
    If keep_ambiguous=True, preserves IUPAC ambiguous bases (RYSWKMBDHVN).
    """
    if keep_ambiguous:
        return re.sub(r"[^ACGTRYSWKMBDHVN]", "", seq.upper())
    else:
        return re.sub(r"[^ATCG]", "", seq.upper())


def validate_dna_sequence(
    sequence: str,
    allow_empty: bool = False,
    allow_ambiguous: bool = False
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate + clean a DNA sequence. Returns (is_valid, cleaned_seq, message_or_None).
    """
    if not sequence and not allow_empty:
        return False, "", "Sequence cannot be empty"

    valid_chars = "ATCG" + ("RYSWKMBDHVN" if allow_ambiguous else "")
    cleaned = "".join(c for c in sequence.upper() if c in valid_chars)
    if not cleaned and sequence:
        return False, "", "Sequence contains no valid DNA characters"
    if len(cleaned) < len(sequence.replace(" ", "")):
        removed = len(sequence.replace(" ", "")) - len(cleaned)
        return True, cleaned, f"Removed {removed} invalid character{'s' if removed != 1 else ''}"
    return True, cleaned, None


# ────────────────────────────────────────────────────────────────────────────────
# 3. Reverse Complement
# ────────────────────────────────────────────────────────────────────────────────
def reverse_complement(seq: str) -> str:
    """
    Compute reverse complement of a DNA string (handles A/T/C/G/N, any case).
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


# ────────────────────────────────────────────────────────────────────────────────
# 4. GC Content
# ────────────────────────────────────────────────────────────────────────────────
def calculate_gc(seq: str) -> float:
    """
    Return GC percentage of a DNA sequence. If seq is empty → 0.0.
    """
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0


# ────────────────────────────────────────────────────────────────────────────────
# 5. Nearest-Neighbor Tₘ Calculation (Average of Breslauer / SantaLucia / Sugimoto)
# ────────────────────────────────────────────────────────────────────────────────
def calculate_tm_consensus(
    sequence: str,
    primer_conc: float = 500e-9,
    na_conc: float = 50e-3
) -> Optional[float]:
    """
    Consensus Tₘ = average of three nearest-neighbor models (Breslauer, SantaLucia, Sugimoto).
    Returns Tₘ in °C, or None if invalid/too short.
    For very short primers (<8 nt), uses a simple approximation: 2°C per A/T, 4°C per G/C minus 7.
    """

    seq = sequence.upper().replace(" ", "")
    if not seq or any(b not in "ATCG" for b in seq):
        return None

    # If < 8 nt → simple Marmur-Doty rule
    if len(seq) < 8:
        a = seq.count("A"); t = seq.count("T")
        g = seq.count("G"); c = seq.count("C")
        return float(2 * (a + t) + 4 * (g + c) - 7)

    R = 1.987  # cal/(mol·K)

    # parameter dictionaries
    breslauer_params = {
        "AA": (-9.1, -24.0), "TT": (-9.1, -24.0),
        "AT": (-8.6, -23.9), "TA": (-6.0, -16.9),
        "CA": (-5.8, -12.9), "TG": (-5.8, -12.9),
        "GT": (-6.5, -17.3), "AC": (-6.5, -17.3),
        "CT": (-7.8, -20.8), "AG": (-7.8, -20.8),
        "GA": (-5.6, -13.5), "TC": (-5.6, -13.5),
        "CG": (-11.9, -27.8), "GC": (-11.1, -26.7),
        "GG": (-11.0, -26.6), "CC": (-11.0, -26.6)
    }

    santalucia_params = {
        "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
        "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
        "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
        "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
        "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
        "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
        "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
        "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)
    }

    sugimoto_params = {
        "AA": (-8.0, -21.9), "TT": (-8.0, -21.9),
        "AT": (-6.6, -17.0), "TA": (-7.2, -19.9),
        "CA": (-8.2, -22.2), "TG": (-8.2, -22.2),
        "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
        "CT": (-7.6, -20.4), "AG": (-7.6, -20.4),
        "GA": (-8.3, -22.2), "TC": (-8.3, -22.2),
        "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
        "GG": (-8.0, -19.9), "CC": (-8.0, -19.9)
    }

    def nn_tm(seq_str: str, params: Dict[str, Tuple[float, float]]) -> float:
        dH = 0.0
        dS = 0.0
        for i in range(len(seq_str) - 1):
            pair = seq_str[i : i + 2]
            if pair in params:
                h, s = params[pair]
                dH += h
                dS += s
            else:
                # fallback for unknown di-nucleotides (very rare)
                dH += -8.0
                dS += -22.0
        # End penalty
        dS += -10.8
        c_eff = primer_conc / 4.0
        tm_kelvin = (dH * 1000.0) / (dS + R * math.log(c_eff))
        tm_c = tm_kelvin - 273.15
        salt_corr = 16.6 * math.log10(na_conc)
        return tm_c + salt_corr

    tm_bres = nn_tm(seq, breslauer_params)
    tm_sant = nn_tm(seq, santalucia_params)
    tm_sugi = nn_tm(seq, sugimoto_params)

    return round((tm_bres + tm_sant + tm_sugi) / 3.0, 2)


# ────────────────────────────────────────────────────────────────────────────────
# 6. Full-Length Primer Design (Automatic Mode)
# ────────────────────────────────────────────────────────────────────────────────
def design_full_length_primers(
    sequence: str,
    min_len: int = 18,
    max_len: int = 30,
    desired_tm: float = 60.0,
    primer_conc: float = 500e-9,
    na_conc: float = 50e-3
) -> Tuple[str, str, float, int, int, float, float]:
    """
    Exhaustively try all forward lengths in [min_len..max_len] at the 5'-end,
    and all reverse lengths in [min_len..max_len] at the 3'-end. Return
      (best_fwd, best_rev, best_score, best_fwd_len, best_rev_len, tm_fwd, tm_rev).
    """
    best_score = float("inf")
    best_fwd = ""
    best_rev = ""
    best_fwd_len = min_len
    best_rev_len = min_len
    best_tm_fwd = 0.0
    best_tm_rev = 0.0

    seq_len = len(sequence)
    if seq_len < (min_len + min_len + 10):
        raise ValueError("Sequence too short for full-length primer design.")

    def score_pair(fwd: str, rev: str) -> Tuple[float, float, float]:
        tm_f = calculate_tm_consensus(fwd, primer_conc, na_conc)
        tm_r = calculate_tm_consensus(rev, primer_conc, na_conc)
        if tm_f is None or tm_r is None:
            return float("inf"), 0.0, 0.0

        tm_diff_penalty = abs(tm_f - tm_r)
        tm_dev_penalty = abs(tm_f - desired_tm) + abs(tm_r - desired_tm)

        def self_dimer(pr: str) -> int:
            rc = reverse_complement(pr)
            score = 0
            for L in range(3, min(7, len(pr) + 1)):
                if pr[-L:] == rc[:L]:
                    score = L
            return score

        dimer_penalty = self_dimer(fwd) + self_dimer(rev)

        def gc_pen(pr: str) -> float:
            gcp = calculate_gc(pr)
            if 40.0 <= gcp <= 60.0:
                return 0.0
            return (40.0 - gcp) if (gcp < 40.0) else (gcp - 60.0)

        gc_penalty = gc_pen(fwd) + gc_pen(rev)

        total = tm_diff_penalty * 1.5 + tm_dev_penalty + dimer_penalty * 2.0 + gc_penalty
        return total, tm_f, tm_r

    for f_len in range(min_len, max_len + 1):
        fwd_cand = sequence[:f_len]
        for r_len in range(min_len, max_len + 1):
            rev_cand = reverse_complement(sequence[-r_len:])
            sc, tm_f, tm_r = score_pair(fwd_cand, rev_cand)
            if sc < best_score:
                best_score = sc
                best_fwd = fwd_cand
                best_rev = rev_cand
                best_fwd_len = f_len
                best_rev_len = r_len
                best_tm_fwd = tm_f
                best_tm_rev = tm_r

    return (
        best_fwd,
        best_rev,
        best_score,
        best_fwd_len,
        best_rev_len,
        best_tm_fwd,
        best_tm_rev,
    )


# ────────────────────────────────────────────────────────────────────────────────
# 7. “Local” Automatic Primer Finder (Used if you want a faster sliding-window 
#     search only within the first/last ~50 nt instead of scanning lengths 18–30)
# ────────────────────────────────────────────────────────────────────────────────
def design_automatic_primers(
    sequence: str,
    fwd_length: int,
    rev_length: int,
    primer_conc: float = 500e-9,
    na_conc: float = 50e-3
) -> Tuple[Tuple[str, str], Tuple[float, float], float, int]:
    """
    Slide a window of size fwd_length in the first ~50 nt, and rev_length in the last ~50 nt,
    score them, pick best. Returns:
      ((best_fwd, best_rev), (tm_fwd, tm_rev), best_score, amplicon_size).
    """

    def score_pair(fwd: str, rev: str) -> Tuple[float, float, float]:
        tm_f = calculate_tm_consensus(fwd, primer_conc, na_conc)
        tm_r = calculate_tm_consensus(rev, primer_conc, na_conc)
        if tm_f is None or tm_r is None:
            return float("inf"), 0.0, 0.0

        tm_diff_penalty = abs(tm_f - tm_r)
        tm_dev_penalty = abs(tm_f - 60.0) + abs(tm_r - 60.0)

        def self_dimer(pr: str) -> int:
            rc = reverse_complement(pr)
            score = 0
            for L in range(3, min(7, len(pr) + 1)):
                if pr[-L:] == rc[:L]:
                    score = L
            return score

        dimer_penalty = self_dimer(fwd) + self_dimer(rev)

        def gc_pen(pr: str) -> float:
            gcp = calculate_gc(pr)
            if 40.0 <= gcp <= 60.0:
                return 0.0
            return (40.0 - gcp) if (gcp < 40.0) else (gcp - 60.0)

        gc_penalty = gc_pen(fwd) + gc_pen(rev)

        total = tm_diff_penalty * 1.5 + tm_dev_penalty + dimer_penalty * 2.0 + gc_penalty
        return total, tm_f, tm_r

    seq_len = len(sequence)
    if seq_len < (fwd_length + rev_length + 10):
        raise ValueError("Sequence too short for primer design.")

    best_score = float("inf")
    best_pair: Tuple[str, str] = ("", "")
    best_tms: Tuple[float, float] = (0.0, 0.0)
    best_amp = 0

    max_rev_window = min(50, seq_len - rev_length)
    for i in range(0, min(50, seq_len - fwd_length + 1)):
        f_cand = sequence[i : i + fwd_length]
        for j in range(max(seq_len - rev_length - max_rev_window, 0), seq_len - rev_length + 1):
            r_cand = reverse_complement(sequence[j : j + rev_length])
            sc, tm_f, tm_r = score_pair(f_cand, r_cand)
            if sc < best_score:
                best_score = sc
                best_pair = (f_cand, r_cand)
                best_tms = (tm_f, tm_r)
                best_amp = (j + rev_length) - i

    return best_pair, best_tms, best_score, best_amp


# ────────────────────────────────────────────────────────────────────────────────
# 8. Primer Design (Cloning)
# ────────────────────────────────────────────────────────────────────────────────
def design_cloning_primers(
    forward_seq: str,
    reverse_seq: str,
    fwd_enzyme: str,
    rev_enzyme: str,
    primer_conc_nM: float = 500,
    custom_prefix: str = "TGCATC"
) -> Tuple[str, str, int, int, Optional[float], Optional[float]]:
    """
    Design cloning primers that sandwich your “inserts” between:
      [custom_prefix] + [fwd_enzyme_site overhang] + [forward insert]
      [custom_prefix] + [rev_enzyme_site overhang] + [reverse-complement(reverse insert)]
    Returns:
      (fwd_primer, rev_primer, fwd_len, rev_len, tm_fwd, tm_rev)
    """
    # Get overhangs from SSD_RESTRICTION_ENZYMES
    over_fwd = SSD_RESTRICTION_ENZYMES.get(fwd_enzyme, {}).get("cut_forward", "")
    over_rev = SSD_RESTRICTION_ENZYMES.get(rev_enzyme, {}).get("cut_forward", "")

    # NdeI special: if forward insert starts with ATG, drop it (its promoter covers it)
    if fwd_enzyme == "NdeI" and forward_seq.upper().startswith("ATG"):
        adj_fwd = forward_seq[3:]
    else:
        adj_fwd = forward_seq

    primer_fwd = custom_prefix + over_fwd + adj_fwd
    primer_rev = custom_prefix + over_rev + reverse_complement(reverse_seq)

    fwd_len = len(primer_fwd)
    rev_len = len(primer_rev)
    tm_fwd = calculate_tm_consensus(primer_fwd, primer_conc_nM * 1e-9)
    tm_rev = calculate_tm_consensus(primer_rev, primer_conc_nM * 1e-9)
    return primer_fwd, primer_rev, fwd_len, rev_len, tm_fwd, tm_rev


# ────────────────────────────────────────────────────────────────────────────────
# 9. (Optional) ORF Finding / Translation / Codon Optimization / etc.
#    – These functions are unchanged from the original, in case other tabs need them.
# ────────────────────────────────────────────────────────────────────────────────

# (a) Translation / Reverse Translation / CODON_USAGE_TABLES / genetic code
GENETIC_CODE: Dict[str, str] = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}

CODON_USAGE_TABLES: Dict[str, Dict[str, List[str]]] = {
    "E. coli BL21": {
        'A': ['GCG', 'GCC', 'GCA', 'GCT'], 'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'],
        'N': ['AAC', 'AAT'], 'D': ['GAT', 'GAC'], 'C': ['TGC', 'TGT'],
        'Q': ['CAG', 'CAA'], 'E': ['GAA', 'GAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'], 'I': ['ATT', 'ATC', 'ATA'],
        'L': ['CTG', 'TTA', 'TTG', 'CTC', 'CTT', 'CTA'], 'K': ['AAA', 'AAG'],
        'M': ['ATG'], 'F': ['TTT', 'TTC'], 'P': ['CCG', 'CCA', 'CCT', 'CCC'],
        'S': ['AGC', 'TCT', 'TCC', 'AGT', 'TCG', 'TCA'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'W': ['TGG'],
        'Y': ['TAT', 'TAC'], 'V': ['GTG', 'GTA', 'GTT', 'GTC'], '*': ['TAA', 'TGA', 'TAG']
    },
    "S. cerevisiae": {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'R': ['AGA', 'AGG', 'CGT', 'CGA', 'CGC', 'CGG'],
        'N': ['AAC', 'AAT'], 'D': ['GAT', 'GAC'], 'C': ['TGT', 'TGC'],
        'Q': ['CAA', 'CAG'], 'E': ['GAA', 'GAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'], 'I': ['ATT', 'ATC', 'ATA'],
        'L': ['TTG', 'CTT', 'TTA', 'CTG', 'CTA', 'CTC'], 'K': ['AAG', 'AAA'],
        'M': ['ATG'], 'F': ['TTT', 'TTC'], 'P': ['CCA', 'CCT', 'CCC', 'CCG'],
        'S': ['TCT', 'TCC', 'TCA', 'AGT', 'TCG', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'W': ['TGG'],
        'Y': ['TAT', 'TAC'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], '*': ['TAA', 'TAG', 'TGA']
    },
    "P. pastoris": {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'R': ['AGA', 'CGT', 'AGG', 'CGA', 'CGC', 'CGG'],
        'N': ['AAC', 'AAT'], 'D': ['GAC', 'GAT'], 'C': ['TGT', 'TGC'],
        'Q': ['CAA', 'CAG'], 'E': ['GAA', 'GAG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'], 'I': ['ATT', 'ATC', 'ATA'],
        'L': ['TTG', 'CTG', 'TTA', 'CTC', 'CTT', 'CTA'], 'K': ['AAG', 'AAA'],
        'M': ['ATG'], 'F': ['TTC', 'TTT'], 'P': ['CCA', 'CCT', 'CCC', 'CCG'],
        'S': ['TCC', 'TCT', 'AGT', 'TCA', 'AGC', 'TCG'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'W': ['TGG'],
        'Y': ['TAC', 'TAT'], 'V': ['GTT', 'GTC', 'GTG', 'GTA'], '*': ['TAA', 'TAG', 'TGA']
    },
    "H. sapiens": {
        'A': ['GCC', 'GCT', 'GCA', 'GCG'], 'R': ['AGG', 'AGA', 'CGG', 'CGC', 'CGA', 'CGT'],
        'N': ['AAC', 'AAT'], 'D': ['GAC', 'GAT'], 'C': ['TGC', 'TGT'],
        'Q': ['CAG', 'CAA'], 'E': ['GAG', 'GAA'], 'G': ['GGC', 'GGG', 'GGA', 'GGT'],
        'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATT', 'ATA'],
        'L': ['CTG', 'CTC', 'TTG', 'CTT', 'TTA', 'CTA'], 'K': ['AAG', 'AAA'],
        'M': ['ATG'], 'F': ['TTC', 'TTT'], 'P': ['CCC', 'CCT', 'CCA', 'CCG'],
        'S': ['AGC', 'TCC', 'TCT', 'AGT', 'TCA', 'TCG'],
        'T': ['ACC', 'ACT', 'ACA', 'ACG'], 'W': ['TGG'],
        'Y': ['TAC', 'TAT'], 'V': ['GTG', 'GTC', 'GTT', 'GTA'], '*': ['TGA', 'TAA', 'TAG']
    },
    "CHO cells": {
        'A': ['GCC', 'GCT', 'GCA', 'GCG'], 'R': ['CGG', 'AGG', 'AGA', 'CGC', 'CGA', 'CGT'],
        'N': ['AAC', 'AAT'], 'D': ['GAC', 'GAT'], 'C': ['TGC', 'TGT'],
        'Q': ['CAG', 'CAA'], 'E': ['GAG', 'GAA'], 'G': ['GGC', 'GGG', 'GGT', 'GGA'],
        'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATT', 'ATA'],
        'L': ['CTG', 'CTC', 'TTG', 'CTT', 'TTA', 'CTA'], 'K': ['AAG', 'AAA'],
        'M': ['ATG'], 'F': ['TTC', 'TTT'], 'P': ['CCC', 'CCT', 'CCA', 'CCG'],
        'S': ['AGC', 'TCC', 'TCT', 'AGT', 'TCA', 'TCG'],
        'T': ['ACC', 'ACT', 'ACA', 'ACG'], 'W': ['TGG'],
        'Y': ['TAC', 'TAT'], 'V': ['GTG', 'GTC', 'GTT', 'GTA'], '*': ['TGA', 'TAA', 'TAG']
    }
}

def translate_sequence(nuc_seq: str, frame: int = 0, find_start: bool = True) -> str:
    """
    Translate DNA → protein (stop at '*' if find_start=True). Uses GENETIC_CODE.
    """
    seq = clean_dna_sequence(nuc_seq)
    if not seq:
        return ""
    if frame not in (0, 1, 2):
        logger.warning(f"Invalid frame {frame}; using frame 0.")
        frame = 0

    start = frame
    if find_start:
        pos = seq.find("ATG", frame)
        if pos != -1:
            start = pos
        else:
            logger.info("No start codon found; translating from frame start.")

    end = len(seq) - ((len(seq) - start) % 3)
    protein = ""
    for i in range(start, end, 3):
        codon = seq[i : i + 3]
        aa = GENETIC_CODE.get(codon, "X")
        protein += aa
        if find_start and aa == "*" and i > start:
            break
    return protein


def reverse_translate_to_dna(prot: str, target_organism: str = "E. coli BL21") -> str:
    """
    Reverse-translate protein → DNA using the most preferred codon from CODON_USAGE_TABLES.
    """
    if target_organism not in CODON_USAGE_TABLES:
        target_organism = "E. coli BL21"
    codon_table = CODON_USAGE_TABLES[target_organism]

    dna = []
    for aa in prot.upper():
        if aa not in codon_table:
            dna.append("NNN")
        else:
            dna.append(codon_table[aa][0])
    return "".join(dna)


# ────────────────────────────────────────────────────────────────────────────────
# 10. (Optional) ORF-Finding / Display
# ────────────────────────────────────────────────────────────────────────────────
STOP_CODONS = ["TAA", "TAG", "TGA"]

def find_orfs(seq: str) -> List[Tuple[int, int, int]]:
    """
    Find all ORFs in DNA. Returns a list of (start, end, frame), 0-based, end is exclusive.
    """
    cleaned = clean_dna_sequence(seq)
    orfs: List[Tuple[int, int, int]] = []

    for frame in range(3):
        i = frame
        while i < len(cleaned) - 2:
            if cleaned[i : i + 3] == "ATG":
                start = i
                for j in range(i + 3, len(cleaned) - 2, 3):
                    if cleaned[j : j + 3] in STOP_CODONS:
                        orfs.append((start, j + 3, frame))
                        i = j + 3
                        break
                else:
                    i += 3
            else:
                i += 1

    return orfs


def display_orfs(seq: str) -> str:
    """
    Return a human-readable list of ORFs with their protein translations.
    """
    orfs = find_orfs(seq)
    if not orfs:
        return "No ORFs found."
    lines = []
    for idx, (start, end, frame) in enumerate(orfs, 1):
        prot = translate_sequence(seq[start:end], 0, False)
        lines.append(f"ORF {idx} (Frame {frame}): {start+1}-{end} bp, Protein: {prot}")
    return "\n".join(lines)


# ────────────────────────────────────────────────────────────────────────────────
# 11. Advanced Codon Optimization (unchanged from original)
# ────────────────────────────────────────────────────────────────────────────────
def advanced_codon_optimization(
    sequence: str,
    target_organism: str = "E. coli BL21",
    optimization_parameters: Optional[Dict] = None,
    is_protein: bool = False
) -> Dict:
    """
    Perform advanced codon optimization preserving logic from G-Synth_2025_5_0.py.
    Returns a dictionary with keys: original_sequence, optimized_sequence, codon_changes, ...
    If an error occurs, returns a dict containing "error".
    """
    if optimization_parameters is None:
        optimization_parameters = {
            'gc_target': (30, 70),
            'avoid_sites': [],
            'avoid_repeats': True,
            'harmonize_usage': True
        }
    results: Dict = {
        "original_sequence": sequence,
        "target_organism": target_organism,
        "is_protein_input": is_protein,
        "optimized_sequence": "",
        "codon_changes": 0,
        "total_codons": 0,
        "gc_before": 0.0,
        "gc_after": 0.0,
        "avoided_sites": [],
        "verification": False
    }
    try:
        # 1) If protein input → reverse translate to DNA
        if is_protein:
            prot_seq = "".join([c for c in sequence.upper() if c.isalpha()])
            dna_seq = reverse_translate_to_dna(prot_seq, target_organism)
            working_seq = dna_seq
            results["total_codons"] = len(prot_seq)
        else:
            cleaned = clean_dna_sequence(sequence)
            working_seq = cleaned
            results["total_codons"] = len(cleaned) // 3

        results["gc_before"] = calculate_gc(working_seq)

        # 2) Optimize codon by codon
        if target_organism not in CODON_USAGE_TABLES:
            target_organism = "E. coli BL21"
        codon_table = CODON_USAGE_TABLES[target_organism]
        optimized = []
        changes = 0

        for i in range(0, len(working_seq), 3):
            codon = working_seq[i : i + 3]
            if len(codon) < 3:
                optimized.append(codon)
                continue
            aa = GENETIC_CODE.get(codon, None)
            if aa is None:
                optimized.append(codon)
                continue

            # 2a) If start codon (i == 0 and AA == 'M'), fix to ATG
            if i == 0 and aa == 'M':
                pick = 'ATG'
                optimized.append(pick)
                if codon != pick:
                    changes += 1
                continue

            # 2b) If stop codon → pick preferred stop from codon_table['*'][0]
            if aa == '*':
                pref_stop = codon_table['*'][0]
                optimized.append(pref_stop)
                if codon != pref_stop:
                    changes += 1
                continue

            # 2c) For non-start, non-stop: get candidate list
            if aa in codon_table:
                candidates = codon_table[aa].copy()
            else:
                optimized.append(codon)
                continue

            # Harmonize usage: if original codon is in the list, keep it
            if optimization_parameters.get('harmonize_usage', True) and codon in candidates:
                pick = codon
            else:
                pick = candidates[0]

            # GC-target adjustment
            gc_min, gc_max = optimization_parameters.get('gc_target', (30, 70))
            current_gc = calculate_gc("".join(optimized))
            if current_gc < gc_min:
                # prefer highest-GC codon
                candidates.sort(key=lambda c: (c.count('G') + c.count('C')), reverse=True)
                pick = candidates[0]
            elif current_gc > gc_max:
                # prefer lowest-GC codon
                candidates.sort(key=lambda c: (c.count('G') + c.count('C')))
                pick = candidates[0]

            # Avoid restriction sites
            avoid_sites = optimization_parameters.get('avoid_sites', [])
            if avoid_sites:
                res_seqs = []
                for site in avoid_sites:
                    rec = SSD_RESTRICTION_ENZYMES.get(site, {}).get("recognition")
                    if rec:
                        res_seqs.append(rec)
                safe_cands = []
                for cand in candidates:
                    context = "".join(optimized[-5:]) + cand + working_seq[i+3 : i+8]
                    conflict = any(rs in context for rs in res_seqs)
                    if not conflict:
                        safe_cands.append(cand)
                if safe_cands:
                    pick = safe_cands[0]
                    # record avoided
                    for rs in res_seqs:
                        if rs in context:
                            results["avoided_sites"].append(rs)

            # Avoid repeats
            if optimization_parameters.get('avoid_repeats', True):
                final_cands = []
                for cand in candidates:
                    ctx = "".join(optimized[-5:]) + cand
                    has_repeat = False
                    for L in range(6, 12):
                        if len(ctx) >= 2 * L:
                            for j in range(len(ctx) - L + 1):
                                fragment = ctx[j : j + L]
                                if ctx.count(fragment) > 1:
                                    has_repeat = True
                                    break
                            if has_repeat:
                                break
                    if not has_repeat:
                        final_cands.append(cand)
                if final_cands:
                    pick = final_cands[0]

            optimized.append(pick)
            if pick != codon:
                changes += 1

        out_seq = "".join(optimized)
        results["optimized_sequence"] = out_seq
        results["codon_changes"] = changes
        results["gc_after"] = calculate_gc(out_seq)

        # Verification: both original & optimized should translate to same protein
        if is_protein:
            opt_pro = translate_sequence(out_seq, 0, False)
            results["verification"] = (opt_pro.replace("*", "") == sequence.replace("*", ""))
        else:
            orig_pro = translate_sequence(working_seq, 0, False)
            opt_pro = translate_sequence(out_seq, 0, False)
            results["verification"] = (orig_pro == opt_pro)

        return results

    except Exception as e:
        logger.error(f"Error in codon optimization: {e}")
        results["error"] = str(e)
        results["optimized_sequence"] = sequence
        return results


# ────────────────────────────────────────────────────────────────────────────────
# 12. Fragmentation for Extended Synthesis (unchanged)
# ────────────────────────────────────────────────────────────────────────────────
ENZYME_PAIRS: Dict[str, Dict[str, str]] = {
    "NdeI / XhoI":  {"forward_overhang": "TA",   "reverse_overhang": "TCGA"},
    "NdeI / EcoRI": {"forward_overhang": "TA",   "reverse_overhang": "AATT"},
    "BamHI / EcoRI":{"forward_overhang": "GATC", "reverse_overhang": "AATT"},
    "BamHI / XhoI":{"forward_overhang": "GATC", "reverse_overhang": "TCGA"},
    "SalI / XbaI":  {"forward_overhang": "TCGAC","reverse_overhang": "TCTAG"}
}

CLEAVAGE_SITES: Dict[str, str] = {
    "Thrombin":    "CTGGTGCCGCGTGGTTCT",
    "TEV":         "GAAAACCTGTATTTTCAGGGC",
    "Factor Xa":   "ATCGAAGGTCGT",
    "PreScission": "CTGGAAGTGCTGTTCCAGGGCCCA",
    "Enterokinase":"GATGACGATGACAAG",
    "SUMO":        "CTGCAGGACTCAGAGG",
    "HRV 3C":      "CTGGAAGTTCTGTTCCAGGGGCCC"
}

def fragment_extended_sequence(
    seq: str,
    frag_size: int,
    enzyme_pair: str,
    cleavage_site: str,
    internal_overlap: int = 15
) -> Tuple[List[Dict], str]:
    """
    BREAK a long sequence into overlapping fragments with sticky ends.
    """
    cleavage_seq = CLEAVAGE_SITES.get(cleavage_site, cleavage_site)
    seq_clean = clean_dna_sequence(seq)
    if len(seq_clean) < frag_size:
        raise ValueError("Sequence is shorter than fragment size.")

    eff_size = frag_size - internal_overlap
    parts: List[str] = []
    for i in range(0, len(seq_clean), eff_size):
        end = min(i + frag_size, len(seq_clean))
        parts.append(seq_clean[i:end])

    assembly = []
    for idx, frag in enumerate(parts):
        if idx == 0:
            fwd = ENZYME_PAIRS[enzyme_pair]["forward_overhang"] + cleavage_seq + frag
            next_ov = parts[idx + 1][:internal_overlap] if idx + 1 < len(parts) else ""
            rev = reverse_complement(frag) + reverse_complement(next_ov)
            typ = "First"
        elif idx == len(parts) - 1:
            prev_ov = parts[idx - 1][-internal_overlap:] if idx > 0 else ""
            fwd = prev_ov + frag
            rev = reverse_complement(frag) + cleavage_seq + ENZYME_PAIRS[enzyme_pair]["reverse_overhang"]
            typ = "Last"
        else:
            prev_ov = parts[idx - 1][-internal_overlap:]
            next_ov = parts[idx + 1][:internal_overlap]
            fwd = prev_ov + frag
            rev = reverse_complement(frag) + reverse_complement(next_ov)
            typ = "Internal"

        assembly.append({
            "fragment": idx + 1,
            "sequence": frag,
            "forward": fwd,
            "reverse": rev,
            "type": typ,
            "length": len(frag)
        })

    reassembled = parts[0]
    for part in parts[1:]:
        reassembled += part[internal_overlap:]
    return assembly, reassembled


# ────────────────────────────────────────────────────────────────────────────────
# 13. Hybridization / Alignment Utilities (unchanged)
# ────────────────────────────────────────────────────────────────────────────────
def is_complement(base1: str, base2: str) -> bool:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return comp.get(base1.upper(), "") == base2.upper()


def optimal_alignment(forward: str, reverse_comp: str, max_shift: Optional[int] = None) -> Tuple[int, int]:
    """
    Find best alignment shift between forward and reverse complement, matching A–T/G–C.
    Returns (best_shift, best_score).
    """
    best_shift = 0
    best_score = 0

    try:
        import numpy as np
        fwd_arr = np.frombuffer(forward.encode(), dtype=np.uint8)
        rev_arr = np.frombuffer(reverse_comp.encode(), dtype=np.uint8)

        match_mat = np.zeros((256, 256), dtype=bool)
        for b1, b2 in [('A','T'),('T','A'),('G','C'),('C','G'),('N','N')]:
            match_mat[ord(b1), ord(b2)] = True

        shifts = range(-len(reverse_comp) + 1, len(forward))
        for shift in shifts:
            o_s = max(0, shift)
            o_e = min(len(forward), shift + len(reverse_comp))
            if o_e <= o_s:
                continue
            f_section = fwd_arr[o_s:o_e]
            r_idx = np.arange(o_s - shift, o_e - shift)
            r_section = rev_arr[r_idx]
            score = int(np.sum(match_mat[f_section, r_section]))
            if score > best_score:
                best_score = score
                best_shift = shift
    except ImportError:
        # Fallback: naive
        for shift in range(-len(reverse_comp) + 1, len(forward)):
            score = 0
            for i in range(max(0, shift), min(len(forward), shift + len(reverse_comp))):
                j = i - shift
                if 0 <= j < len(reverse_comp) and is_complement(forward[i], reverse_comp[j]):
                    score += 1
            if score > best_score:
                best_score = score
                best_shift = shift

    return best_shift, best_score


def validate_sequence_complementarity(seq1: str, seq2: str) -> Tuple[bool, List[Tuple[int, str, str]]]:
    """
    Check if seq2 is the reverse complement of seq1 (same length).
    Returns (is_complementary, list_of_mismatches).
    """
    if len(seq1) != len(seq2):
        return False, []

    mismatches = []
    for i, (b1, b2) in enumerate(zip(seq1.upper(), seq2.upper())):
        if not is_complement(b1, b2):
            mismatches.append((i, b1, b2))
    return (len(mismatches) == 0), mismatches


# ────────────────────────────────────────────────────────────────────────────────
# 14. One-letter / Three-letter AA Conversion (unchanged)
# ────────────────────────────────────────────────────────────────────────────────
THREE_TO_ONE: Dict[str, str] = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Stop': '*'
}

ONE_TO_THREE: Dict[str, str] = {v: k for k, v in THREE_TO_ONE.items()}


def convert_to_three_letter(prot: str) -> str:
    """
    Convert one-letter protein sequence → three-letter codes (space separated).
    """
    return " ".join(ONE_TO_THREE.get(aa, aa) for aa in prot)


def convert_three_to_one(prot: str) -> str:
    """
    Convert three-letter codes (space separated) → one-letter string.
    """
    res = ""
    for tok in prot.split():
        res += THREE_TO_ONE.get(tok, "")
    return res


# ────────────────────────────────────────────────────────────────────────────────
# 15. Optional: Biopython Compatibility (if installed)
# ────────────────────────────────────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    USING_BIOPYTHON = True
except ImportError:
    USING_BIOPYTHON = False
    logger.warning("Biopython not found; some features will be limited.")


def load_sequence_from_fasta(contents: str) -> str:
    """
    Given FASTA text (">header\nATGC..."), return cleaned concatenated DNA sequence.
    """
    lines = contents.strip().splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return clean_dna_sequence(seq)


def seq_to_biopython_record(seq: str, record_id: str = "seq1"):
    """
    If Biopython is installed, return a SeqRecord for downstream visualization.
    """
    if not USING_BIOPYTHON:
        raise ImportError("Biopython not installed.")
    from Bio.SeqRecord import SeqRecord
    return SeqRecord(Seq(seq), id=record_id)

# ────────────────────────────────────────────────────────────────────────────────
# End of utils/bio_utils.py
# ────────────────────────────────────────────────────────────────────────────────
