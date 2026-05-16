"""Pairwise & multiple sequence alignment."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import streamlit as st
from gsynth_ui.state import init_state
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.alignment import align_global, align_local, multiple_alignment
from gsynth_core.io import parse_fasta

st.set_page_config(page_title="Alignment — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("📏 Sequence Alignment")

tab_pair, tab_msa = st.tabs(["🔗 Pairwise", "🧮 Multiple (MSA)"])

with tab_pair:
    col1, col2 = st.columns(2)
    with col1:
        seq_a = st.text_area("Sequence A", height=120, key="pa_a")
    with col2:
        seq_b = st.text_area("Sequence B", height=120, key="pa_b")
    mode = st.radio("Mode", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"], horizontal=True)
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        match = st.number_input("Match", value=2.0)
    with col2:
        mismatch = st.number_input("Mismatch", value=-1.0)
    with col3:
        gap_open = st.number_input("Gap open", value=-10.0)
    with col4:
        gap_extend = st.number_input("Gap extend", value=-0.5)
    if st.button("Align", type="primary", key="pa_run"):
        a = "".join(seq_a.split()).upper()
        b = "".join(seq_b.split()).upper()
        if not a or not b:
            st.error("Both sequences required.")
        else:
            fn = align_global if mode.startswith("Global") else align_local
            aln = fn(a, b, match=match, mismatch=mismatch,
                     gap_open=gap_open, gap_extend=gap_extend)
            c1, c2, c3 = st.columns(3)
            with c1: st.markdown(stat("Score", f"{aln.score:.1f}"), unsafe_allow_html=True)
            with c2: st.markdown(stat("Identity", f"{aln.identity_percent:.1f}%"), unsafe_allow_html=True)
            with c3: st.markdown(stat("Mode", aln.mode), unsafe_allow_html=True)
            st.markdown(f"<div class='gs-mono'>{aln.pretty()}</div>", unsafe_allow_html=True)

with tab_msa:
    st.markdown("**Paste FASTA with 2+ sequences, or upload a file.**")
    msa_text = st.text_area("FASTA", height=200, key="msa_text",
                             placeholder=">a\nATGAAA…\n>b\nATGAAG…\n")
    upload = st.file_uploader("Or upload", type=["fasta", "fa", "fna"], key="msa_up")
    if upload is not None:
        msa_text = upload.read().decode("utf-8", errors="replace")
    method = st.selectbox("Method", ["auto", "mafft", "muscle", "clustalo", "star"])
    if st.button("Align all", type="primary", key="msa_run"):
        if not msa_text.strip():
            st.error("Need FASTA input.")
        else:
            try:
                records = parse_fasta(msa_text)
                if len(records) < 2:
                    st.error("Need at least 2 sequences.")
                else:
                    pairs = [(r.name or f"seq{i+1}", r.sequence) for i, r in enumerate(records)]
                    msa = multiple_alignment(pairs, method=method)
                    c1, c2, c3 = st.columns(3)
                    with c1: st.markdown(stat("Sequences", str(len(msa.aligned))), unsafe_allow_html=True)
                    with c2: st.markdown(stat("Columns", str(msa.length)), unsafe_allow_html=True)
                    with c3: st.markdown(stat("Method", msa.method), unsafe_allow_html=True)
                    st.markdown(f"<div class='gs-mono'>{msa.pretty()}</div>", unsafe_allow_html=True)
            except Exception as e:
                st.error(f"Alignment failed: {e}")
