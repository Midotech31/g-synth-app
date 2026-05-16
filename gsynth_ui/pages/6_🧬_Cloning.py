"""Cloning page — Gibson, Golden Gate, restriction."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
import pandas as pd
from gsynth_ui.state import init_state
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.cloning import (
    design_gibson_overlaps, simulate_gibson_assembly,
    design_golden_gate, simulate_golden_gate_assembly,
    simulate_restriction_cloning,
)
from gsynth_core.restriction import list_enzymes

st.set_page_config(page_title="Cloning — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("🧬 In-silico Cloning")

tab_gibson, tab_gg, tab_restr = st.tabs(["🪡 Gibson Assembly", "🧩 Golden Gate", "✂ Restriction Cloning"])

with tab_gibson:
    st.markdown("**Enter 2-10 fragments. Adjacent fragments must already share an overlap region.**")
    n_frags = st.number_input("Number of fragments", 2, 10, 3, key="gib_n")
    fragments: list[tuple[str, str]] = []
    for i in range(n_frags):
        col1, col2 = st.columns([1, 4])
        with col1:
            name = st.text_input(f"Name #{i+1}", value=f"f{i+1}", key=f"gib_name_{i}")
        with col2:
            seq = st.text_area(f"Sequence #{i+1}", height=80, key=f"gib_seq_{i}",
                                placeholder="ATGCATGC…")
        fragments.append((name, "".join(seq.split()).upper()))
    target_tm = st.slider("Target overlap Tm (°C)", 40.0, 65.0, 50.0, key="gib_tm")
    circular = st.checkbox("Circular assembly (plasmid)", True, key="gib_circ")
    if st.button("Design + Simulate", type="primary", key="gib_run"):
        if not all(seq for _, seq in fragments):
            st.error("All fragments need a sequence.")
        else:
            try:
                designed = design_gibson_overlaps(
                    fragments, target_tm=target_tm, circular=circular,
                )
                df = pd.DataFrame([{
                    "name": f.name, "length": len(f.sequence),
                    "left overlap": f.left_overlap, "left Tm": f"{f.left_tm:.1f}",
                    "right overlap": f.right_overlap, "right Tm": f"{f.right_tm:.1f}",
                } for f in designed])
                st.dataframe(df, use_container_width=True, hide_index=True)
                assembled = simulate_gibson_assembly(designed, circular=circular)
                st.markdown(f"### Assembled product ({len(assembled)} nt)")
                st.markdown(f"<div class='gs-mono'>{assembled}</div>", unsafe_allow_html=True)
            except Exception as e:
                st.error(f"Failed: {e}")

with tab_gg:
    st.markdown("**Enter parts as inserts only — flanking BsaI/BsmBI sites are added automatically.**")
    n_parts = st.number_input("Number of parts", 2, 8, 4, key="gg_n")
    enz = st.selectbox("Type-IIS enzyme", ["BsaI", "BsmBI"], key="gg_enz")
    parts: list[tuple[str, str]] = []
    for i in range(n_parts):
        col1, col2 = st.columns([1, 4])
        with col1:
            name = st.text_input(f"Name #{i+1}", value=f"p{i+1}", key=f"gg_name_{i}")
        with col2:
            ins = st.text_area(f"Insert #{i+1}", height=70, key=f"gg_ins_{i}")
        parts.append((name, "".join(ins.split()).upper()))
    if st.button("Design + Simulate", type="primary", key="gg_run"):
        if not all(ins for _, ins in parts):
            st.error("All parts need a sequence.")
        else:
            try:
                gg = design_golden_gate(parts, enzyme=enz)
                df = pd.DataFrame([{
                    "name": p.name, "left fusion": p.left_fusion,
                    "right fusion": p.right_fusion, "insert len": len(p.insert),
                    "full len (with enzyme sites)": len(p.full_sequence),
                } for p in gg])
                st.dataframe(df, use_container_width=True, hide_index=True)
                assembled = simulate_golden_gate_assembly(gg, circular=True)
                st.markdown(f"### Assembled circular product ({len(assembled)} nt)")
                st.markdown(f"<div class='gs-mono'>{assembled}</div>", unsafe_allow_html=True)
                if "GGTCTC" in assembled or "GAGACC" in assembled:
                    st.error("BsaI site leaked into the product!")
                else:
                    st.success(f"✅ No {enz} site in the final product — assembly is irreversible.")
            except Exception as e:
                st.error(f"Failed: {e}")

with tab_restr:
    st.markdown("**Classical double-digest + ligation.**")
    enzymes = list_enzymes("curated")
    col1, col2 = st.columns(2)
    with col1:
        ea = st.selectbox("Enzyme A", enzymes,
                          index=enzymes.index("EcoRI") if "EcoRI" in enzymes else 0, key="rc_ea")
    with col2:
        eb = st.selectbox("Enzyme B", enzymes,
                          index=enzymes.index("BamHI") if "BamHI" in enzymes else 0, key="rc_eb")
    vector = st.text_area("Vector sequence (circular)", height=150, key="rc_vec")
    insert = st.text_area("Insert sequence (PCR product, with flanking enzyme sites)",
                          height=120, key="rc_ins")
    if st.button("Simulate cloning", type="primary", key="rc_run"):
        vec = "".join(vector.split()).upper()
        ins = "".join(insert.split()).upper()
        if not vec or not ins:
            st.error("Need both vector and insert.")
        else:
            try:
                clone = simulate_restriction_cloning(vec, ins, enzyme_a=ea, enzyme_b=eb)
                if clone is None:
                    st.warning("No compatible clone found — overhangs do not align.")
                else:
                    c1, c2, c3 = st.columns(3)
                    with c1: st.markdown(stat("Length", f"{len(clone.sequence)} nt"), unsafe_allow_html=True)
                    with c2: st.markdown(stat("Backbone", f"{clone.vector_fragment.length} nt"), unsafe_allow_html=True)
                    with c3: st.markdown(stat("Insert", f"{clone.insert_fragment.length} nt"), unsafe_allow_html=True)
                    st.markdown(f"**Directional:** {'✅' if clone.directional else '⚠ no — self-ligation possible'}")
                    st.markdown("### Final clone")
                    st.markdown(f"<div class='gs-mono'>{clone.sequence}</div>", unsafe_allow_html=True)
            except Exception as e:
                st.error(f"Failed: {e}")
