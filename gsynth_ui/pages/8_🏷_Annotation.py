"""Plasmid auto-annotation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
import pandas as pd
from gsynth_ui.state import init_state, get_active_sequence
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.annotation import annotate_features, KNOWN_FEATURES
from gsynth_core.io import write_genbank

st.set_page_config(page_title="Annotation — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("🏷 Plasmid Auto-annotation")

active = get_active_sequence()
if active is None:
    st.info("👈 Load a sequence first.")
    st.stop()
seq = active.sequence
st.caption(f"Active: **{active.name}** ({len(seq):,} nt)")

min_id = st.slider("Minimum identity for fuzzy match", 0.5, 1.0, 0.85, 0.05)
append = st.checkbox("Add detected features to the record", True)

if st.button("Scan", type="primary"):
    hits = annotate_features(active, min_identity=min_id, append_to_record=append)
    c1, c2 = st.columns(2)
    with c1: st.markdown(stat("Hits", str(len(hits))), unsafe_allow_html=True)
    with c2: st.markdown(stat("DB size", str(len(KNOWN_FEATURES))), unsafe_allow_html=True)
    if hits:
        df = pd.DataFrame([{
            "feature": h.name, "type": h.type,
            "start": h.start, "end": h.end, "strand": h.strand,
            "identity": f"{h.identity:.0%}",
        } for h in hits])
        st.dataframe(df, use_container_width=True, hide_index=True)

        if append:
            gb = write_genbank(active)
            st.download_button("Download annotated GenBank", gb,
                                file_name=f"{active.name}.gb", mime="text/plain")
    else:
        st.info("No known features detected. Try lowering the identity threshold.")

with st.expander("📚 Feature database"):
    db = pd.DataFrame([{"name": k, "type": v[0], "length": len(v[1])}
                        for k, v in KNOWN_FEATURES.items()])
    st.dataframe(db, use_container_width=True, hide_index=True)
