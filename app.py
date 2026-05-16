"""G-Synth 3.0 — Streamlit application entry point.

Run with: `streamlit run app.py`
"""
from __future__ import annotations

import sys
from pathlib import Path

# Ensure repo root is on path
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import streamlit as st

from gsynth_ui.state import init_state, get_sequences, get_active_sequence, set_active, add_sequence, remove_sequence
from gsynth_ui.theme import apply_theme, stat
from gsynth_core import __version__
from gsynth_core.io import parse_fasta, parse_genbank
from gsynth_core.io.records import SeqRecord
from gsynth_core.sequence import gc_content


st.set_page_config(
    page_title="G-Synth 3.0",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

init_state()
apply_theme()


# ── Sidebar ────────────────────────────────────────────────────────────────
_LOGO_PATH = _HERE / "assets" / "logo.png"

with st.sidebar:
    if _LOGO_PATH.exists():
        st.image(str(_LOGO_PATH), use_container_width=True)
    st.markdown(f"## 🧬 G-Synth\n*v{__version__}*")

    st.markdown("### 📂 Open sequence")
    upload = st.file_uploader(
        "Upload FASTA or GenBank",
        type=["fasta", "fa", "fna", "gb", "gbk"],
        label_visibility="collapsed",
    )
    if upload is not None:
        text = upload.read().decode("utf-8", errors="replace")
        try:
            if upload.name.lower().endswith((".gb", ".gbk")):
                records = parse_genbank(text)
            else:
                records = parse_fasta(text)
            for r in records:
                if not r.name:
                    r.name = upload.name.rsplit(".", 1)[0]
                add_sequence(r)
            st.success(f"Loaded {len(records)} record(s)")
        except Exception as e:
            st.error(f"Failed to parse: {e}")

    st.markdown("**Paste sequence**")
    raw = st.text_area("Sequence", key="paste_seq", height=120,
                       placeholder="ATGAAACGC… (FASTA, plain DNA, or protein)",
                       label_visibility="collapsed")
    paste_name = st.text_input("Name", key="paste_name", value="pasted", label_visibility="collapsed",
                                placeholder="Sequence name")
    if st.button("Add pasted", use_container_width=True) and raw.strip():
        try:
            if raw.lstrip().startswith(">"):
                for r in parse_fasta(raw):
                    add_sequence(r)
            else:
                rec = SeqRecord(sequence="".join(raw.split()).upper(), name=paste_name or "pasted")
                add_sequence(rec)
            st.rerun()
        except Exception as e:
            st.error(f"Could not parse: {e}")

    st.markdown("### 📑 Sequences")
    seqs = get_sequences()
    if not seqs:
        st.caption("_No sequences loaded._")
    else:
        active = st.session_state.get("active", "")
        names = list(seqs.keys())
        sel = st.radio("Active", names, index=names.index(active) if active in names else 0,
                       label_visibility="collapsed")
        if sel != active:
            set_active(sel)
            st.rerun()
        if st.button("🗑 Remove active", use_container_width=True):
            remove_sequence(active)
            st.rerun()

    st.divider()
    st.caption("Use the page navigator above (Streamlit sidebar) to switch modules.")


# ── Main page: dashboard ──────────────────────────────────────────────────
st.title("🧬 G-Synth 3.0 — Dashboard")

active = get_active_sequence()
if active is None:
    st.info(
        "👈 Load a sequence from the sidebar to begin.\n\n"
        "**Try one of these sample workflows:**\n"
        "- Load a FASTA/GenBank file or paste a sequence,\n"
        "- Then pick a module in the **pages/** navigator (top-left).\n"
    )
    st.markdown("---")
    st.markdown("### Available modules")
    st.markdown(
        "- **Sequence Analysis** — translation, ORF finder, GC profile\n"
        "- **Codon Optimization** — multi-objective CAI + GC + restriction-site avoidance\n"
        "- **Primer Design** — PCR primer pairs with optional restriction-site overhangs\n"
        "- **Restriction Sites** — 80-enzyme search + digestion simulator\n"
        "- **CRISPR Designer** — Doench 2016 + CFD off-target\n"
        "- **Cloning** — Gibson, Golden Gate, classical restriction\n"
        "- **Plasmid Annotation** — auto-detect promoters, ORIs, resistance markers, tags\n"
        "- **Pairwise Alignment** — Needleman-Wunsch + Smith-Waterman\n"
    )
else:
    seq = active.sequence
    gc = gc_content(seq)
    st.markdown(f"### {active.name}")
    if active.description:
        st.caption(active.description)
    st.markdown(
        f"<div class='gs-card'>"
        f"{stat('Length', f'{len(seq):,} nt')}"
        f"{stat('GC content', f'{gc:.1f}%')}"
        f"{stat('Topology', active.topology)}"
        f"{stat('Features', f'{len(active.features)}')}"
        f"</div>",
        unsafe_allow_html=True,
    )
    st.markdown("**Sequence preview**")
    preview = seq[:200] + (" …" if len(seq) > 200 else "")
    st.markdown(f"<div class='gs-mono'>{preview}</div>", unsafe_allow_html=True)
    if active.features:
        st.markdown("### Features")
        feat_rows = [
            {"type": f.type, "label": f.name, "start": f.location.start,
             "end": f.location.end, "strand": f.location.strand}
            for f in active.features
        ]
        st.dataframe(feat_rows, use_container_width=True, hide_index=True)


st.markdown(
    "<div class='gs-footer'>G-Synth 3.0 — Dr. Mohamed Merzoug, ESSBO · MIT License</div>",
    unsafe_allow_html=True,
)
