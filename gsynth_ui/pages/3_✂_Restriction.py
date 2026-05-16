"""Restriction site search and digest simulation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
import pandas as pd
from gsynth_ui.state import init_state, get_active_sequence
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.restriction import (
    find_sites, simulate_double_digestion, list_enzymes, suggest_compatible_ends,
)

st.set_page_config(page_title="Restriction — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("✂ Restriction Sites & Digestion")

active = get_active_sequence()
if active is None:
    st.info("👈 Load a sequence first.")
    st.stop()
seq = active.sequence
st.caption(f"Active: **{active.name}** ({len(seq):,} nt, {active.topology})")

tab_search, tab_digest, tab_compat = st.tabs(["🔍 Find sites", "🧪 Digest", "🔗 Compatible ends"])

with tab_search:
    enzyme_filter = st.multiselect("Enzymes (empty = all 80 curated)",
                                    list_enzymes("curated"))
    sites = find_sites(seq, enzymes=enzyme_filter or None)
    st.markdown(f"**{len(sites)} site{'s' if len(sites) != 1 else ''} found**")
    if sites:
        df = pd.DataFrame([{
            "enzyme": s.enzyme.name, "site": s.enzyme.site, "position": s.position,
            "strand": s.strand, "overhang": s.enzyme.overhang_sequence or "(blunt)",
            "type": s.enzyme.overhang_type.value, "cut@fwd": s.cut_position_fwd,
            "cut@rev": s.cut_position_rev,
        } for s in sites])
        st.dataframe(df, use_container_width=True, hide_index=True)
        # Frequency-by-enzyme summary
        freq = pd.DataFrame(df.groupby("enzyme").size().reset_index(name="count").sort_values("count", ascending=False))
        st.markdown("**Frequencies**")
        st.dataframe(freq, use_container_width=True, hide_index=True)

with tab_digest:
    enzymes = list_enzymes("curated")
    col1, col2 = st.columns(2)
    with col1:
        ea = st.selectbox("Enzyme A", enzymes, index=enzymes.index("EcoRI") if "EcoRI" in enzymes else 0)
    with col2:
        eb = st.selectbox("Enzyme B", enzymes, index=enzymes.index("BamHI") if "BamHI" in enzymes else 0)
    circular = st.checkbox("Treat as circular plasmid", active.topology == "circular")
    if st.button("Simulate digest", type="primary"):
        frags = simulate_double_digestion(seq, ea, eb, circular=circular)
        c1, c2 = st.columns(2)
        with c1: st.markdown(stat("Fragments", str(len(frags))), unsafe_allow_html=True)
        with c2: st.markdown(stat("Total bp", f"{sum(f.length for f in frags):,}"), unsafe_allow_html=True)
        df = pd.DataFrame([{
            "#": i + 1, "length": f.length,
            "left overhang": f.left_overhang or "(blunt)",
            "left type": f.left_overhang_type.value,
            "right overhang": f.right_overhang or "(blunt)",
            "right type": f.right_overhang_type.value,
            "start@template": f.start_on_template, "end@template": f.end_on_template,
        } for i, f in enumerate(frags)])
        st.dataframe(df, use_container_width=True, hide_index=True)

with tab_compat:
    enzymes = list_enzymes("curated")
    target = st.selectbox("Find compatible ends for", enzymes,
                          index=enzymes.index("BamHI") if "BamHI" in enzymes else 0)
    compat = suggest_compatible_ends(target)
    if compat:
        st.markdown(f"**{target}** is compatible with: {', '.join(compat)}")
    else:
        st.info(f"No commercially-curated enzyme produces an overhang compatible with {target}.")
