"""Primer design page."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
import pandas as pd
from gsynth_ui.state import init_state, get_active_sequence
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.primers import design_primer_pair, design_cloning_primers, analyze_primer, PrimerParams
from gsynth_core.restriction import list_enzymes

st.set_page_config(page_title="Primers — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("🧪 Primer Design")

active = get_active_sequence()
if active is None:
    st.info("👈 Load a sequence first.")
    st.stop()
seq = active.sequence
st.caption(f"Active: **{active.name}** ({len(seq):,} nt)")

with st.form("primer_form"):
    col1, col2, col3 = st.columns(3)
    with col1:
        opt_tm = st.number_input("Optimal Tm (°C)", 50.0, 75.0, 60.0)
        max_dtm = st.number_input("Max ΔTm (°C)", 0.5, 10.0, 3.0)
    with col2:
        min_len = st.number_input("Min length", 15, 40, 18)
        max_len = st.number_input("Max length", 18, 60, 30)
    with col3:
        cloning = st.checkbox("Add cloning sites (restriction)", False)
        if cloning:
            enzymes = list_enzymes("curated")
            le = st.selectbox("5' enzyme (forward)", enzymes,
                               index=enzymes.index("EcoRI") if "EcoRI" in enzymes else 0)
            re_enz = st.selectbox("3' enzyme (reverse)", enzymes,
                                   index=enzymes.index("BamHI") if "BamHI" in enzymes else 0)
    submit = st.form_submit_button("🚀 Design primers", type="primary", use_container_width=True)

if submit:
    params = PrimerParams(min_len=min_len, max_len=max_len, optimal_tm=opt_tm, max_tm_diff=max_dtm)
    if cloning:
        pair = design_cloning_primers(seq, left_enzyme=le, right_enzyme=re_enz, params=params)
    else:
        pair = design_primer_pair(seq, params=params)
    if pair is None:
        st.error("No primer pair satisfies the constraints. Try relaxing Tm range or length.")
    else:
        c1, c2, c3 = st.columns(3)
        with c1: st.markdown(stat("Product", f"{pair.product_size} bp"), unsafe_allow_html=True)
        with c2: st.markdown(stat("ΔTm", f"{pair.tm_difference:.2f} °C"), unsafe_allow_html=True)
        with c3: st.markdown(stat("Score", f"{pair.score:.1f}"), unsafe_allow_html=True)

        st.markdown("### Forward primer")
        an_fwd = analyze_primer(pair.forward.sequence)
        st.markdown(f"<div class='gs-mono'>{pair.forward.sequence}</div>", unsafe_allow_html=True)
        st.markdown(f"Length={pair.forward.length}, Tm={pair.forward.tm:.1f}°C, "
                    f"GC={pair.forward.gc:.0f}%, hairpin risk={an_fwd.hairpin.risk_level}")

        st.markdown("### Reverse primer")
        an_rev = analyze_primer(pair.reverse.sequence)
        st.markdown(f"<div class='gs-mono'>{pair.reverse.sequence}</div>", unsafe_allow_html=True)
        st.markdown(f"Length={pair.reverse.length}, Tm={pair.reverse.tm:.1f}°C, "
                    f"GC={pair.reverse.gc:.0f}%, hairpin risk={an_rev.hairpin.risk_level}")
