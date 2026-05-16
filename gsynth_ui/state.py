"""Centralized session-state management.

Streamlit re-runs the script on every interaction. To keep sequences and
results across reruns, every module reads/writes via this module's helpers.
"""
from __future__ import annotations

import streamlit as st

from gsynth_core.io.records import SeqRecord


_KEYS = {
    "sequences": dict,   # name → SeqRecord
    "active": str,       # name of the currently selected sequence
    "results": dict,     # module name → last result object
    "ui_theme": str,
}


def init_state() -> None:
    """Initialise all expected keys on first app run."""
    for key, kind in _KEYS.items():
        if key not in st.session_state:
            st.session_state[key] = kind() if kind is not str else ""
    if "ui_theme" not in st.session_state or not st.session_state["ui_theme"]:
        st.session_state["ui_theme"] = "default"


def get_sequences() -> dict[str, SeqRecord]:
    return st.session_state.get("sequences", {})


def get_active_sequence() -> SeqRecord | None:
    seqs = get_sequences()
    name = st.session_state.get("active", "")
    return seqs.get(name)


def add_sequence(record: SeqRecord, *, make_active: bool = True) -> None:
    seqs = get_sequences()
    base_name = record.name or "untitled"
    name = base_name
    i = 2
    while name in seqs:
        name = f"{base_name}_{i}"
        i += 1
    record.name = name
    seqs[name] = record
    st.session_state["sequences"] = seqs
    if make_active:
        st.session_state["active"] = name


def remove_sequence(name: str) -> None:
    seqs = get_sequences()
    if name in seqs:
        del seqs[name]
        if st.session_state.get("active") == name:
            st.session_state["active"] = next(iter(seqs), "")
    st.session_state["sequences"] = seqs


def set_active(name: str) -> None:
    if name in get_sequences():
        st.session_state["active"] = name


def store_result(module: str, result) -> None:
    res = st.session_state.get("results", {})
    res[module] = result
    st.session_state["results"] = res


def get_result(module: str):
    return st.session_state.get("results", {}).get(module)
