"""CRISPR guide RNA design."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
import pandas as pd
from gsynth_ui.state import init_state, get_active_sequence
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.crispr import design_guides, CasType

st.set_page_config(page_title="CRISPR — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("✨ CRISPR Guide Designer")
st.caption("Doench 2016 on-target + CFD off-target scoring")

active = get_active_sequence()
if active is None:
    st.info("👈 Load a sequence first.")
    st.stop()
seq = active.sequence
st.caption(f"Active: **{active.name}** ({len(seq):,} nt)")

with st.form("crispr_form"):
    col1, col2, col3 = st.columns(3)
    with col1:
        cas = st.selectbox("Cas variant",
                            ["SpCas9", "SpCas9-NG", "SaCas9", "Cas12a", "Cas13"])
    with col2:
        min_on = st.slider("Min on-target score", 0.0, 1.0, 0.4, 0.05)
    with col3:
        min_gc = st.slider("Min GC %", 0, 100, 30)
        max_gc = st.slider("Max GC %", 0, 100, 70)
    both = st.checkbox("Both strands", True)
    use_ref = st.checkbox("Use this sequence as off-target reference", False)
    top_n = st.number_input("Top N", 1, 50, 10)
    submit = st.form_submit_button("🔍 Find guides", type="primary", use_container_width=True)

if submit:
    cas_enum = CasType(cas)
    guides = design_guides(
        seq, cas=cas_enum, both_strands=both,
        min_gc=min_gc, max_gc=max_gc, min_on_target=min_on,
        genome_reference=seq if use_ref else None,
        top_n=top_n,
    )
    st.markdown(f"**{len(guides)} guides found**")
    if guides:
        rows = []
        for g in guides:
            rows.append({
                "protospacer": g.protospacer, "PAM": g.pam,
                "position": g.start, "strand": g.strand,
                "GC %": f"{g.gc_percent:.0f}",
                "on-target": f"{g.on_target:.2f}",
                "specificity": f"{g.specificity:.1f}" if g.specificity is not None else "—",
                "top OT CFD": f"{g.top_off_target_cfd:.3f}" if g.top_off_target_cfd is not None else "—",
                "OT count": g.off_target_count if g.off_target_count is not None else "—",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
