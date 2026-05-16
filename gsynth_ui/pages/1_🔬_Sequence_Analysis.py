"""Sequence analysis page — translation, ORFs, GC profile."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
import pandas as pd
from gsynth_ui.state import init_state, get_active_sequence
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.sequence import six_frame_translate, find_orfs, gc_content, reverse_complement

st.set_page_config(page_title="Sequence Analysis — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("🔬 Sequence Analysis")

active = get_active_sequence()
if active is None:
    st.info("👈 Load a sequence first.")
    st.stop()

seq = active.sequence
st.caption(f"Active: **{active.name}** ({len(seq):,} nt)")

c1, c2, c3 = st.columns(3)
with c1: st.markdown(stat("Length", f"{len(seq):,} nt"), unsafe_allow_html=True)
with c2: st.markdown(stat("GC %", f"{gc_content(seq):.1f}"), unsafe_allow_html=True)
with c3: st.markdown(stat("Reverse complement", f"{len(reverse_complement(seq)):,} nt"), unsafe_allow_html=True)

tab_orf, tab_translate, tab_rc = st.tabs(["🧬 ORF finder", "📜 6-frame translation", "↔ Reverse complement"])

with tab_orf:
    min_len = st.number_input("Min ORF length (aa)", 10, 1000, 50)
    both = st.checkbox("Both strands", True)
    alt = st.checkbox("Allow alternative start codons (GTG, TTG) — prokaryotic", False)
    orfs = find_orfs(seq, min_length_aa=min_len, both_strands=both, allow_alt_starts=alt)
    st.markdown(f"**{len(orfs)} ORF{'s' if len(orfs) != 1 else ''} found**")
    if orfs:
        df = pd.DataFrame([{
            "start": o.start, "end": o.end, "strand": o.strand, "frame": o.frame,
            "length (aa)": o.length_aa, "start codon": o.start_codon,
            "stop codon": o.stop_codon or "—",
            "protein (preview)": o.protein[:30] + ("…" if len(o.protein) > 30 else ""),
        } for o in orfs])
        st.dataframe(df, use_container_width=True, hide_index=True)

with tab_translate:
    frames = six_frame_translate(seq)
    for f in ("+1", "+2", "+3", "-1", "-2", "-3"):
        st.markdown(f"**Frame {f}**")
        st.markdown(f"<div class='gs-mono'>{frames[f]}</div>", unsafe_allow_html=True)

with tab_rc:
    rc = reverse_complement(seq)
    st.markdown(f"<div class='gs-mono'>{rc}</div>", unsafe_allow_html=True)
    st.download_button("Download as FASTA",
                       f">{active.name}_rc\n{rc}\n",
                       file_name=f"{active.name}_rc.fasta", mime="text/plain")
