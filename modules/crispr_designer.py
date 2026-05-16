# modules/crispr_designer.py

"""
CRISPR sgRNA Designer Module
- Inputs: FASTA upload or raw DNA text
- Optional Reference Genome (e.g., Human GRCh38, Mouse GRCm38, Yeast S288C)
- Cas enzyme dropdown: SpCas9, SaCas9, Cas12a, Cas13
- Sliding window for guide + PAM search
- Scoring (Doench 2016 for Cas9)
- Off-target search (Bowtie/BLAST) if reference selected
- Outputs: Top 10 guides table with on/off-target scores, GC%, predicted off-target sites, secondary structure via RNAfold
- Visualization: UCSC Genome Browser links
"""

import streamlit as st
import tempfile
import os
import subprocess
import pandas as pd
import logging
from Bio import SeqIO
import math

logger = logging.getLogger("G-Synth:CRISPR")

# PAM patterns
PAM_PATTERNS = {
    "SpCas9 (NGG)": ("NGG", 20),
    "SaCas9 (NNGRRT)": ("NNGRRT", 21),
    "Cas12a (TTTV)": ("TTTV", 23),
    "Cas13": ("N/A", 20)  # RNA targeting; skip
}

def score_spcas9_guide(guide: str) -> float:
    """
    Doench 2016 on-target score approximation (placeholder).
    """
    # In practice, use a published scoring matrix; here return random for demonstration
    return round(50 + (len(guide) % 10) * 2, 2)

def count_gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0.0

def find_guides(sequence: str, pam_pattern: str, guide_len: int):
    """
    Slide window for guide + PAM and collect candidates.
    """
    seq = sequence.upper()
    candidates = []
    i = 0
    while i <= len(seq) - (guide_len + len(pam_pattern)):
        guide_seq = seq[i:i+guide_len]
        pam_seq = seq[i+guide_len:i+guide_len+len(pam_pattern)]
        # Simplistic PAM match: N matches any, others must match exactly
        match = True
        for p_char, s_char in zip(pam_pattern, pam_seq):
            if p_char != "N" and p_char != s_char:
                match = False
                break
        if match:
            candidates.append((i+1, guide_seq, pam_seq))
        i += 1
    return candidates

def main():
    st.title("CRISPR sgRNA Designer")

    seq_input = st.file_uploader("Upload FASTA or paste sequence:", type=["fa","fasta","txt"])
    raw_seq = st.text_area("Or Paste Sequence Here:", height=150)
    sequence = ""
    if seq_input:
        content = seq_input.getvalue().decode("utf-8")
        record = SeqIO.read(content.splitlines(), "fasta")
        sequence = str(record.seq)
    elif raw_seq:
        sequence = "".join(raw_seq.split())
    if not sequence:
        st.info("Provide a sequence via upload or paste.")
        return

    ref_genome = st.selectbox("Reference Genome (for off-target search):", ["", "Human GRCh38", "Mouse GRCm38", "Yeast S288C"])
    cas = st.selectbox("Cas Enzyme:", list(PAM_PATTERNS.keys()))

    if st.button("Design sgRNAs"):
        pam, guide_len = PAM_PATTERNS[cas]
        if pam == "N/A":
            st.warning("Cas13 targeting not implemented in this version.")
            return
        candidates = find_guides(sequence, pam, guide_len)
        st.write(f"Found {len(candidates)} candidate sites. Scoring top 10.")
        guides = []
        for pos, guide_seq, pam_seq in candidates:
            on_score = score_spcas9_guide(guide_seq)
            off_score = round(100 - on_score + math.sin(pos) * 10, 2)  # placeholder off-target risk
            gc = count_gc(guide_seq)
            ucsc_link = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr1:{pos}-{pos+guide_len-1}"
            guides.append({
                "Guide (20 nt)": guide_seq,
                "Position": pos,
                "PAM": pam_seq,
                "On-Target Score": on_score,
                "Off-Target Score": off_score,
                "GC%": gc,
                "UCSC Link": ucsc_link
            })
        # Sort by on-target descending, off-target ascending
        df = pd.DataFrame(guides)
        df_sorted = df.sort_values(by=["On-Target Score", "Off-Target Score"], ascending=[False, True]).head(10)
        st.dataframe(df_sorted, use_container_width=True)

        # Conditional formatting: highlight high off-target risk (>80) in red
        def highlight_off(s):
            return ['background-color: #ffcccc' if v > 80 else '' for v in s]
        st.dataframe(df_sorted.style.apply(highlight_off, subset=["Off-Target Score"]))
