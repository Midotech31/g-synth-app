"""Single source of truth for UI styling — no per-module CSS allowed."""
from __future__ import annotations
import streamlit as st


_CSS = """
<style>
  :root {
    --gs-primary: #1e3a8a;
    --gs-accent: #06b6d4;
    --gs-bg-card: #f8fafc;
    --gs-border: #e2e8f0;
    --gs-text: #0f172a;
    --gs-text-muted: #64748b;
    --gs-radius: 8px;
  }
  .gs-card {
    background: var(--gs-bg-card);
    border: 1px solid var(--gs-border);
    border-radius: var(--gs-radius);
    padding: 1rem;
    margin-bottom: 1rem;
  }
  .gs-card h3 {
    margin-top: 0;
    color: var(--gs-primary);
    font-size: 1.05rem;
    font-weight: 600;
  }
  .gs-stat {
    display: inline-block;
    margin-right: 1.5rem;
  }
  .gs-stat .label {
    font-size: 0.8rem;
    color: var(--gs-text-muted);
  }
  .gs-stat .value {
    font-size: 1.3rem;
    font-weight: 600;
    color: var(--gs-text);
  }
  .gs-pill {
    display: inline-block;
    padding: 0.15rem 0.5rem;
    border-radius: 999px;
    background: var(--gs-accent);
    color: white;
    font-size: 0.75rem;
    font-weight: 500;
    margin-right: 0.25rem;
  }
  .gs-mono {
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
    font-size: 0.85rem;
    line-height: 1.4;
    color: var(--gs-text);
    background: #ffffff;
    border: 1px solid var(--gs-border);
    border-radius: var(--gs-radius);
    padding: 0.5rem;
    overflow-x: auto;
    white-space: pre;
  }
  .gs-footer {
    color: var(--gs-text-muted);
    font-size: 0.8rem;
    text-align: center;
    padding: 1rem 0;
  }
</style>
"""


def apply_theme() -> None:
    """Inject CSS once per page (must be called from each page)."""
    st.markdown(_CSS, unsafe_allow_html=True)


def stat(label: str, value: str) -> str:
    return f"""<div class="gs-stat"><div class="label">{label}</div>
<div class="value">{value}</div></div>"""


def card(title: str, body_html: str) -> str:
    return f"""<div class="gs-card"><h3>{title}</h3>{body_html}</div>"""
