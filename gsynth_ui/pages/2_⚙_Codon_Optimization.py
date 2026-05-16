"""Codon optimization page."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import streamlit as st
from gsynth_ui.state import init_state, get_active_sequence, store_result
from gsynth_ui.theme import apply_theme, stat
from gsynth_core.codon import optimize, OptimizationParams
from gsynth_core.sequence import translate, gc_content

st.set_page_config(page_title="Codon Optimization — G-Synth", layout="wide")
init_state()
apply_theme()

st.title("⚙ Codon Optimization")
st.caption("Multi-objective: CAI + GC window + homopolymer + repeat + restriction-site avoidance")

active = get_active_sequence()

with st.form("opt_form"):
    col1, col2 = st.columns([2, 1])
    with col1:
        default_seq = ""
        if active is not None:
            default_seq = active.sequence
        protein_or_dna = st.text_area("Protein OR DNA sequence", default_seq, height=180,
                                       help="Auto-detected: amino-acid letters → protein, ACGT → DNA (will be translated)")
    with col2:
        organism = st.selectbox("Target organism", ["e_coli", "human", "yeast", "pichia", "cho"])
        gc_min = st.number_input("Min GC (%)", 0.0, 100.0, 40.0)
        gc_max = st.number_input("Max GC (%)", 0.0, 100.0, 65.0)
        max_homo = st.number_input("Max homopolymer", 1, 10, 5)
        forbid_str = st.text_input("Forbidden sites (comma-separated)",
                                    "GAATTC,GGATCC,AAGCTT,CTCGAG")
    submit = st.form_submit_button("🚀 Optimize", use_container_width=True, type="primary")

if submit:
    raw = "".join(protein_or_dna.split()).upper()
    if not raw:
        st.error("Empty sequence")
    else:
        looks_like_dna = set(raw) <= set("ACGTN")
        if looks_like_dna:
            try:
                protein = translate(raw, to_stop=False).rstrip("*")
                st.info(f"DNA detected — translated to {len(protein)} aa.")
            except Exception as e:
                st.error(f"Translation failed: {e}")
                st.stop()
            original_dna: str | None = raw
        else:
            protein = raw
            original_dna = None

        forbid = tuple(s.strip().upper() for s in forbid_str.split(",") if s.strip())
        try:
            res = optimize(
                protein,
                params=OptimizationParams(
                    organism=organism, gc_window_min=gc_min, gc_window_max=gc_max,
                    avoid_sites=forbid, max_homopolymer=max_homo,
                ),
                original_dna=original_dna,
            )
        except Exception as e:
            st.error(f"Optimization failed: {e}")
            st.stop()

        store_result("codon_optimize", res)
        st.success("Optimization complete!")

        c1, c2, c3, c4 = st.columns(4)
        with c1: st.markdown(stat("CAI", f"{res.cai_after:.3f}"), unsafe_allow_html=True)
        with c2: st.markdown(stat("GC %", f"{res.gc_after:.1f}"), unsafe_allow_html=True)
        with c3: st.markdown(stat("Length", f"{len(res.optimized_dna)} nt"), unsafe_allow_html=True)
        with c4: st.markdown(stat("Sites OK", "✅" if not res.forbidden_sites_remaining else "❌"),
                              unsafe_allow_html=True)

        if res.cai_before is not None:
            st.markdown(f"**Before optimization:** CAI={res.cai_before:.3f}, "
                        f"GC={res.gc_before:.1f}%, changes={res.changes} nt")

        if res.forbidden_sites_remaining:
            st.warning(f"Could not remove: {', '.join(res.forbidden_sites_remaining)}")
        for w in res.warnings:
            st.warning(w)

        st.markdown("### Optimized DNA")
        st.markdown(f"<div class='gs-mono'>{res.optimized_dna}</div>", unsafe_allow_html=True)
        st.download_button("Download FASTA",
                            f">{organism}_optimized\n{res.optimized_dna}\n",
                            file_name="optimized.fasta", mime="text/plain")
