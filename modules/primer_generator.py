# -*- coding: utf-8 -*-
"""
Enhanced Primer Generator Module with Fixed UI Issues and Proper Highlighting
"""

import re
from typing import Tuple
import streamlit as st

# Import core functions from bio_utils
from utils.bio_utils import (
    reverse_complement,
    calculate_gc,
    calculate_tm_consensus,
    design_full_length_primers,
    design_cloning_primers,
)

# Import restriction enzymes
try:
    from modules.ssd import SSD_RESTRICTION_ENZYMES
except ImportError:
    SSD_RESTRICTION_ENZYMES = {}

def inject_primer_css():
    """Enhanced CSS with fixed highlighting and layout issues"""
    st.markdown("""
    <style>
    /* ═══════════════════════════════════════════════════════════════════════════════
       FIXED PRIMER GENERATOR STYLING WITH PROPER HIGHLIGHTING
    ═══════════════════════════════════════════════════════════════════════════════ */
    
    /* Main header styling */
    .primer-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 2rem;
        border-radius: 16px;
        margin-bottom: 2rem;
        text-align: center;
        box-shadow: 0 8px 32px rgba(102, 126, 234, 0.3);
        border: 3px solid rgba(255, 255, 255, 0.2);
    }
    
    .primer-header h1 {
        margin: 0;
        font-size: 2.5rem;
        font-weight: 800;
        text-shadow: 0 2px 4px rgba(0,0,0,0.3);
        letter-spacing: -1px;
    }
    
    .primer-header p {
        margin: 0.5rem 0 0 0;
        font-size: 1.1rem;
        opacity: 0.9;
        font-weight: 500;
    }
    
    /* Mode selector styling */
    .mode-selector {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
        border: 2px solid #e2e8f0;
        border-radius: 16px;
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        transition: all 0.3s ease;
    }
    
    .mode-selector:hover {
        box-shadow: 0 8px 25px rgba(0,0,0,0.12);
        transform: translateY(-2px);
    }
    
    /* Enhanced input containers */
    .input-container {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
        border: 2px solid #e2e8f0;
        border-radius: 16px;
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        transition: all 0.3s ease;
    }
    
    .input-container:hover {
        border-color: #667eea;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15);
        transform: translateY(-2px);
    }
    
    /* Results container */
    .results-container {
        background: linear-gradient(135deg, #f0fff4 0%, #f7fafc 100%);
        border: 2px solid #10b981;
        border-radius: 16px;
        padding: 2rem;
        margin: 2rem 0;
        box-shadow: 0 8px 25px rgba(16, 185, 129, 0.15);
        position: relative;
        overflow: hidden;
    }
    
    .results-container::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        height: 4px;
        background: linear-gradient(90deg, #10b981, #059669, #047857);
    }
    
    /* FIXED: Primer card styling with proper header containment */
    .primer-card {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
        border: 2px solid #e2e8f0;
        border-radius: 12px;
        padding: 0;
        margin: 1rem 0;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        transition: all 0.3s ease;
        overflow: hidden;
        position: relative;
    }
    
    .primer-card:hover {
        border-color: #667eea;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15);
        transform: translateY(-4px);
    }
    
    /* FIXED: Primer card header (colored band/ruban) */
    .primer-card-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 0.75rem 1.5rem;
        margin: 0;
        font-size: 1.1rem;
        font-weight: 600;
        text-align: center;
        border-bottom: 2px solid rgba(255, 255, 255, 0.2);
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    
    /* FIXED: Primer card content area */
    .primer-card-content {
        padding: 1.5rem;
    }
    
    /* FIXED: Enhanced primer sequence styling with proper highlighting support */
    .primer-sequence {
        font-family: 'Courier New', 'Monaco', 'Menlo', monospace !important;
        font-size: 1rem !important;
        font-weight: bold !important;
        background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%) !important;
        padding: 1rem !important;
        border-radius: 8px !important;
        border-left: 4px solid #667eea !important;
        margin: 1rem 0 !important;
        word-break: break-all !important;
        line-height: 1.8 !important;
        letter-spacing: 1px !important;
        border: 2px solid #e2e8f0 !important;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.06) !important;
        transition: all 0.3s ease !important;
        white-space: pre-wrap !important;
        overflow-wrap: break-word !important;
    }
    
    .primer-sequence:hover {
        background: linear-gradient(135deg, #e2e8f0 0%, #cbd5e1 100%) !important;
        border-color: #667eea !important;
    }
    
    /* FIXED: Proper highlighting for cloning primers - Enhanced visibility */
    .highlight-prefix {
        background: linear-gradient(135deg, #FFD700, #FFC107) !important;
        padding: 3px 6px !important;
        border-radius: 4px !important;
        font-weight: bold !important;
        border: 2px solid #DAA520 !important;
        box-shadow: 0 2px 4px rgba(218, 165, 32, 0.5) !important;
        color: #8B4513 !important;
        text-shadow: 0 1px 2px rgba(255,255,255,0.8) !important;
        margin: 0 1px !important;
        display: inline !important;
    }
    
    .highlight-enzyme {
        background: linear-gradient(135deg, #FFFF00, #FFED4E) !important;
        padding: 3px 6px !important;
        border-radius: 4px !important;
        font-weight: bold !important;
        border: 2px solid #CCCC00 !important;
        box-shadow: 0 2px 4px rgba(204, 204, 0, 0.5) !important;
        color: #B45309 !important;
        text-shadow: 0 1px 2px rgba(255,255,255,0.8) !important;
        margin: 0 1px !important;
        display: inline !important;
    }
    
    /* FIXED: Component breakdown styling */
    .component-info {
        background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 0.75rem 1rem;
        margin: 0.5rem 0;
        font-family: 'Courier New', monospace;
        font-size: 0.9rem;
        border-left: 4px solid #667eea;
        transition: all 0.3s ease;
    }
    
    .component-info:hover {
        background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%);
        transform: translateX(4px);
    }
    
    /* Cloning info panel */
    .cloning-info {
        background: linear-gradient(135deg, #eff6ff 0%, #dbeafe 100%);
        border: 2px solid #3b82f6;
        border-radius: 12px;
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 15px rgba(59, 130, 246, 0.15);
    }
    
    /* Property badges */
    .property-badge {
        display: inline-block;
        background: linear-gradient(135deg, #667eea, #764ba2);
        color: white;
        padding: 0.4rem 0.8rem;
        border-radius: 20px;
        font-size: 0.85rem;
        font-weight: 600;
        margin: 0.2rem;
        box-shadow: 0 2px 8px rgba(102, 126, 234, 0.3);
        border: 1px solid rgba(255, 255, 255, 0.2);
    }
    
    /* Status indicators */
    .status-success {
        background: linear-gradient(135deg, #10b981, #059669);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 20px;
        font-weight: 600;
        display: inline-block;
        margin: 0.5rem 0;
        box-shadow: 0 2px 8px rgba(16, 185, 129, 0.3);
    }
    
    .status-warning {
        background: linear-gradient(135deg, #f59e0b, #d97706);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 20px;
        font-weight: 600;
        display: inline-block;
        margin: 0.5rem 0;
        box-shadow: 0 2px 8px rgba(245, 158, 11, 0.3);
    }
    
    /* Enhanced buttons */
    .stButton > button {
        background: linear-gradient(135deg, #667eea, #764ba2) !important;
        color: white !important;
        border: none !important;
        border-radius: 12px !important;
        padding: 0.75rem 2rem !important;
        font-weight: 600 !important;
        font-size: 1.1rem !important;
        transition: all 0.3s ease !important;
        box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3) !important;
        width: 100% !important;
        border: 2px solid rgba(255, 255, 255, 0.2) !important;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px) !important;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.4) !important;
        background: linear-gradient(135deg, #5a67d8, #667eea) !important;
    }
    
    /* Enhanced input fields */
    .stTextInput > div > div > input,
    .stTextArea > div > div > textarea,
    .stSelectbox > div > div > select,
    .stNumberInput > div > div > input {
        border-radius: 8px !important;
        border: 2px solid #e2e8f0 !important;
        transition: all 0.3s ease !important;
        font-size: 0.95rem !important;
    }
    
    .stTextInput > div > div > input:focus,
    .stTextArea > div > div > textarea:focus,
    .stSelectbox > div > div > select:focus,
    .stNumberInput > div > div > input:focus {
        border-color: #667eea !important;
        box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1) !important;
    }
    
    /* Metrics enhancement */
    .stMetric {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
        padding: 1rem;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
        transition: all 0.3s ease;
    }
    
    .stMetric:hover {
        box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        transform: translateY(-2px);
    }
    
    /* Info panels */
    .info-panel {
        background: linear-gradient(135deg, #eff6ff 0%, #dbeafe 100%);
        border-left: 4px solid #3b82f6;
        padding: 1.5rem;
        border-radius: 12px;
        margin: 1rem 0;
        box-shadow: 0 4px 15px rgba(59, 130, 246, 0.15);
    }
    
    .info-panel h4 {
        color: #1e40af;
        margin: 0 0 1rem 0;
        font-weight: 700;
    }
    
    .info-panel ul {
        margin: 0;
        color: #1e40af;
    }
    
    /* Responsive design */
    @media (max-width: 768px) {
        .primer-header h1 {
            font-size: 2rem;
        }
        .primer-card-content {
            padding: 1rem;
        }
        .primer-sequence {
            font-size: 0.9rem !important;
            letter-spacing: 0.5px !important;
        }
    }
    </style>
    """, unsafe_allow_html=True)

def render_enhanced_header():
    """Render modern header for primer generator"""
    st.markdown("""
    <div class="primer-header">
        <h1>🧬 Primer Generator</h1>
        <p>Advanced PCR & Cloning Primer Design Platform</p>
    </div>
    """, unsafe_allow_html=True)

def render_mode_selector():
    """Enhanced mode selector with interactive buttons"""
    st.markdown('<div class="mode-selector">', unsafe_allow_html=True)
    st.markdown("### 🎯 Select Design Mode")
    
    # Mode descriptions
    mode_descriptions = {
        "Manual": "📝 Direct primer extraction from sequence ends",
        "Automatic": "🤖 AI-optimized primer pair selection", 
        "Cloning": "🔬 Restriction enzyme cloning primers"
    }
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("📝 Manual", key="mode_manual", use_container_width=True):
            st.session_state.primer_mode = "Manual"
    with col2:
        if st.button("🤖 Automatic", key="mode_auto", use_container_width=True):
            st.session_state.primer_mode = "Automatic"
    with col3:
        if st.button("🔬 Cloning", key="mode_cloning", use_container_width=True):
            st.session_state.primer_mode = "Cloning"
    
    # Initialize mode if not set
    if 'primer_mode' not in st.session_state:
        st.session_state.primer_mode = "Manual"
    
    # Display current mode with description
    current_mode = st.session_state.primer_mode
    st.markdown(f"""
    <div style="text-align: center; margin-top: 1rem;">
        <span class="status-success">Active: {current_mode} Mode</span>
        <p style="margin: 0.5rem 0; color: #6b7280; font-style: italic;">
            {mode_descriptions[current_mode]}
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)
    
    return current_mode

def render_sequence_input():
    """Enhanced sequence input with real-time validation"""
    st.markdown('<div class="input-container">', unsafe_allow_html=True)
    st.markdown("### 🧬 DNA Sequence Input")
    
    raw_seq = st.text_area(
        "Enter DNA Sequence (5′→3′)",
        height=140,
        placeholder="ATCGATCGATCG...",
        help="Enter nucleotide sequence using A, T, C, G, N characters only"
    )
    
    # Real-time sequence validation
    if raw_seq:
        sequence = re.sub(r"[^ATCGN]", "", raw_seq.upper())
        invalid_chars = len(raw_seq) - len(sequence)
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Sequence Length", f"{len(sequence)} nt")
        with col2:
            gc_content = calculate_gc(sequence) if sequence else 0
            st.metric("GC Content", f"{gc_content:.1f}%")
        with col3:
            if invalid_chars > 0:
                st.metric("Invalid Chars", invalid_chars, delta=-invalid_chars)
            else:
                st.metric("Status", "✅ Valid")
    else:
        sequence = ""
    
    st.markdown('</div>', unsafe_allow_html=True)
    return sequence

def render_primer_results(results):
    """FIXED: Enhanced results display with proper layout and highlighting"""
    if not results:
        return
        
    st.markdown('<div class="results-container">', unsafe_allow_html=True)
    st.markdown("## 🎯 Primer Design Results")
    
    col1, col2 = st.columns(2, gap="large")
    
    # FIXED: Forward primer card with proper header containment
    with col1:
        st.markdown("""
        <div class="primer-card">
            <div class="primer-card-header">➡️ Forward Primer</div>
            <div class="primer-card-content">
        """, unsafe_allow_html=True)
        
        if results["mode"] == "Cloning":
            highlighted = highlight_cloning_sequence(
                results["fwd"], results["enz_fwd"], results["prefix"]
            )
            st.markdown(f'<div class="primer-sequence">{highlighted}</div>', 
                       unsafe_allow_html=True)
            
            # Add component breakdown for cloning primers
            st.markdown("**Components:**")
            if results["prefix"] in results["fwd"]:
                st.markdown(f'<div class="component-info">🟡 <strong>Stuffer:</strong> <code>{results["prefix"]}</code></div>', 
                           unsafe_allow_html=True)
            if results["enz_fwd"] in SSD_RESTRICTION_ENZYMES:
                recog = SSD_RESTRICTION_ENZYMES[results["enz_fwd"]]["recognition"]
                if recog in results["fwd"]:
                    st.markdown(f'<div class="component-info">🟨 <strong>{results["enz_fwd"]} Site:</strong> <code>{recog}</code></div>', 
                               unsafe_allow_html=True)
        else:
            st.markdown(f'<div class="primer-sequence">{results["fwd"]}</div>', 
                       unsafe_allow_html=True)
        
        # Properties as badges
        st.markdown(f"""
        <div>
            <span class="property-badge">Length: {results['len_f']} nt</span>
            <span class="property-badge">GC: {results['gc_f']:.1f}%</span>
            <span class="property-badge">Tm: {results['tm_f']:.1f}°C</span>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown('</div></div>', unsafe_allow_html=True)
    
    # FIXED: Reverse primer card with proper header containment
    with col2:
        st.markdown("""
        <div class="primer-card">
            <div class="primer-card-header">⬅️ Reverse Primer</div>
            <div class="primer-card-content">
        """, unsafe_allow_html=True)
        
        if results["mode"] == "Cloning":
            highlighted = highlight_cloning_sequence(
                results["rev"], results["enz_rev"], results["prefix"]
            )
            st.markdown(f'<div class="primer-sequence">{highlighted}</div>', 
                       unsafe_allow_html=True)
            
            # Add component breakdown for cloning primers
            st.markdown("**Components:**")
            if results["prefix"] in results["rev"]:
                st.markdown(f'<div class="component-info">🟡 <strong>Stuffer:</strong> <code>{results["prefix"]}</code></div>', 
                           unsafe_allow_html=True)
            if results["enz_rev"] in SSD_RESTRICTION_ENZYMES:
                recog = SSD_RESTRICTION_ENZYMES[results["enz_rev"]]["recognition"]
                if recog in results["rev"]:
                    st.markdown(f'<div class="component-info">🟨 <strong>{results["enz_rev"]} Site:</strong> <code>{recog}</code></div>', 
                               unsafe_allow_html=True)
        else:
            st.markdown(f'<div class="primer-sequence">{results["rev"]}</div>', 
                       unsafe_allow_html=True)
        
        # Properties as badges
        st.markdown(f"""
        <div>
            <span class="property-badge">Length: {results['len_r']} nt</span>
            <span class="property-badge">GC: {results['gc_r']:.1f}%</span>
            <span class="property-badge">Tm: {results['tm_r']:.1f}°C</span>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown('</div></div>', unsafe_allow_html=True)
    
    # Cloning-specific information
    if results["mode"] == "Cloning":
        st.markdown("---")
        st.markdown('<div class="cloning-info">', unsafe_allow_html=True)
        st.markdown("### 🔬 Cloning Information")
        col3, col4, col5 = st.columns(3)
        with col3:
            st.markdown(f"**Forward Enzyme:** {results['enz_fwd']}")
            if results["enz_fwd"] in SSD_RESTRICTION_ENZYMES:
                recog = SSD_RESTRICTION_ENZYMES[results["enz_fwd"]]["recognition"]
                st.markdown(f"Recognition: `{recog}`")
        with col4:
            st.markdown(f"**Reverse Enzyme:** {results['enz_rev']}")
            if results["enz_rev"] in SSD_RESTRICTION_ENZYMES:
                recog = SSD_RESTRICTION_ENZYMES[results["enz_rev"]]["recognition"]
                st.markdown(f"Recognition: `{recog}`")
        with col5:
            st.markdown(f"**Stuffer Prefix:** `{results['prefix']}`")
            st.markdown(f"Length: {len(results['prefix'])} nt")
        st.markdown('</div>', unsafe_allow_html=True)
    
    # General metrics
    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Amplicon Size", f"{results['size']} bp")
    with col2:
        tm_diff = abs(results['tm_f'] - results['tm_r'])
        st.metric("Tm Difference", f"{tm_diff:.1f}°C", delta=f"{tm_diff:.1f}")
    with col3:
        avg_gc = (results['gc_f'] + results['gc_r']) / 2
        st.metric("Average GC", f"{avg_gc:.1f}%")
    
    # Copy to clipboard functionality
    if st.button("📋 Copy Results to Clipboard", key="copy_results"):
        clip_text = (
            f"Forward Primer: {results['fwd']}\n"
            f"Reverse Primer: {results['rev']}\n"
            f"Forward: {results['len_f']} nt, GC {results['gc_f']:.1f}%, Tm {results['tm_f']:.1f}°C\n"
            f"Reverse: {results['len_r']} nt, GC {results['gc_r']:.1f}%, Tm {results['tm_r']:.1f}°C\n"
            f"Amplicon Size: {results['size']} bp"
        )
        st.code(clip_text, language="text")
        st.success("✅ Results copied! Use Ctrl+C (Cmd+C) to copy the text above.")
    
    st.markdown('</div>', unsafe_allow_html=True)

# ──────────────────────────────────────────────────────────────────────────────
# Helper Functions (Complete Implementation - All Business Logic Preserved)
# ──────────────────────────────────────────────────────────────────────────────

def compute_self_dimer_score(primer: str) -> int:
    """Compute self-dimerization score"""
    rc = reverse_complement(primer)
    best = 0
    for length in range(3, min(7, len(primer) + 1)):
        if primer[-length:] == rc[:length]:
            best = length
    return best

def compute_gc_penalty(primer: str, ideal_range: Tuple[float, float] = (40.0, 60.0)) -> float:
    """Return GC penalty"""
    gc = calculate_gc(primer)
    low, high = ideal_range
    if low <= gc <= high:
        return 0.0
    return (low - gc) if (gc < low) else (gc - high)

def score_primer_pair(fwd: str, rev: str, desired_tm: float = 60.0, 
                     primer_conc: float = 500e-9, na_conc: float = 50e-3) -> Tuple[float, float, float]:
    """Score primer pair"""
    tm_fwd = calculate_tm_consensus(fwd, primer_conc, na_conc)
    tm_rev = calculate_tm_consensus(rev, primer_conc, na_conc)

    if tm_fwd is None or tm_rev is None:
        return float("inf"), 0.0, 0.0

    tm_diff_penalty = abs(tm_fwd - tm_rev) * 1.5
    tm_dev_penalty = abs(tm_fwd - desired_tm) + abs(tm_rev - desired_tm)
    dimer_penalty = (compute_self_dimer_score(fwd) + compute_self_dimer_score(rev)) * 2.0
    gc_penalty = compute_gc_penalty(fwd) + compute_gc_penalty(rev)

    total = tm_diff_penalty + tm_dev_penalty + dimer_penalty + gc_penalty
    return total, tm_fwd, tm_rev

def design_automatic_primers(sequence: str, fwd_len: int, rev_len: int, 
                           primer_conc: float = 500e-9) -> Tuple[str, str, float, float, int]:
    """Design automatic primers"""
    best_score = float("inf")
    best_values = None
    fwd_window = rev_window = 50

    if len(sequence) < fwd_len + rev_len + 10:
        raise ValueError("Sequence too short for chosen primer lengths.")

    for i in range(0, min(fwd_window, len(sequence) - fwd_len + 1)):
        fwd_candidate = sequence[i : i + fwd_len]
        for j in range(max(len(sequence) - rev_len - rev_window, 0), len(sequence) - rev_len + 1):
            rev_candidate = reverse_complement(sequence[j : j + rev_len])
            score, tm_fwd, tm_rev = score_primer_pair(fwd_candidate, rev_candidate, primer_conc=primer_conc)
            if score < best_score:
                best_score = score
                amplicon_size = (j + rev_len) - i
                best_values = (fwd_candidate, rev_candidate, tm_fwd, tm_rev, amplicon_size)

    if best_values is None:
        raise RuntimeError("Could not find any valid primer pair.")
    return best_values

def highlight_cloning_sequence(primer: str, enzyme: str, prefix: str = "TGCATC") -> str:
    """
    FIXED: Highlight cloning sequence components with enhanced HTML formatting
    Returns HTML with highlighted stuffer (prefix) and restriction site
    """
    if enzyme not in SSD_RESTRICTION_ENZYMES:
        return primer

    recog = SSD_RESTRICTION_ENZYMES[enzyme]["recognition"]
    seq = primer
    
    # Find the prefix (stuffer) position
    pre_pos = seq.find(prefix)
    if pre_pos == -1:
        return primer
    
    # Find the restriction site position after the prefix
    site_pos = seq.find(recog, pre_pos + len(prefix))
    if site_pos == -1:
        return primer

    # Build the highlighted sequence
    before = seq[:pre_pos]
    between = seq[pre_pos + len(prefix) : site_pos]
    after = seq[site_pos + len(recog) :]
    
    highlighted = (
        f"{before}"
        f"<span class='highlight-prefix'>{prefix}</span>"
        f"{between}"
        f"<span class='highlight-enzyme'>{recog}</span>"
        f"{after}"
    )
    
    return highlighted

def main() -> None:
    """Complete enhanced main function with all fixes and beautiful design"""
    # Inject enhanced CSS
    inject_primer_css()
    
    # Render enhanced header
    render_enhanced_header()
    
    # Enhanced mode selector
    mode = render_mode_selector()
    
    # Enhanced sequence input
    sequence = render_sequence_input()
    
    if not sequence:
        st.info("👆 Please enter a DNA sequence to begin primer design")
        return
    
    # Enhanced parameters section
    st.markdown('<div class="input-container">', unsafe_allow_html=True)
    st.markdown("### ⚙️ Design Parameters")
    
    col1, col2 = st.columns(2)
    with col1:
        primer_conc_nM = st.number_input(
            "Primer Concentration (nM)", 
            min_value=10, max_value=2000, value=500, step=10,
            help="Typical range: 100-1000 nM"
        )
        fwd_len = st.number_input(
            "Forward Primer Length (nt)", 
            min_value=5, max_value=100, value=20, step=1
        )
    with col2:
        st.write("")  # Spacing
        rev_len = st.number_input(
            "Reverse Primer Length (nt)", 
            min_value=5, max_value=100, value=20, step=1
        )
    
    primer_conc = primer_conc_nM * 1e-9
    
    # Mode-specific options
    if mode == "Automatic":
        full_amplicon = st.checkbox(
            "🔍 Full-Length Amplicon Analysis",
            help="Search across entire sequence with variable primer lengths (18-30 nt)"
        )
    else:
        full_amplicon = False
    
    if mode == "Cloning":
        enz_list = list(SSD_RESTRICTION_ENZYMES.keys())
        col3, col4 = st.columns(2)
        with col3:
            enz_fwd = st.selectbox("Forward Enzyme", enz_list, 
                                 index=enz_list.index("NdeI") if "NdeI" in enz_list else 0,
                                 help="Choose restriction enzyme for forward primer")
        with col4:
            enz_rev = st.selectbox("Reverse Enzyme", enz_list,
                                 index=enz_list.index("XhoI") if "XhoI" in enz_list else 0,
                                 help="Choose restriction enzyme for reverse primer")
        
        use_prefix = st.checkbox("🧪 Custom Stuffer Prefix")
        prefix = st.text_input("Stuffer Prefix", value="TGCATC") if use_prefix else "TGCATC"
    else:
        enz_fwd = enz_rev = ""
        prefix = "TGCATC"
    
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Generate button
    if st.button("🚀 Generate Primers", key="generate"):
        with st.spinner("🧬 Designing optimal primers..."):
            try:
                # All existing primer generation logic (unchanged)
                if mode == "Manual":
                    fwd = sequence[:fwd_len]
                    rev = reverse_complement(sequence[-rev_len:])
                    tm_fwd = calculate_tm_consensus(fwd, primer_conc)
                    tm_rev = calculate_tm_consensus(rev, primer_conc)
                    amplicon_size = len(sequence)
                    fwd_len_eff, rev_len_eff = len(fwd), len(rev)

                elif mode == "Automatic":
                    if full_amplicon:
                        (fwd, rev, _score, _lf, _lr, tm_fwd, tm_rev) = design_full_length_primers(
                            sequence, 18, 30, 60.0, primer_conc
                        )
                        amplicon_size = len(sequence)
                        fwd_len_eff, rev_len_eff = len(fwd), len(rev)
                    else:
                        (fwd, rev, tm_fwd, tm_rev, amplicon_size) = design_automatic_primers(
                            sequence, fwd_len, rev_len, primer_conc
                        )
                        fwd_len_eff, rev_len_eff = fwd_len, rev_len

                else:  # Cloning mode
                    insert_fwd = sequence[:fwd_len]
                    insert_rev = sequence[-rev_len:]
                    (fwd, rev, fwd_len_eff, rev_len_eff, tm_fwd, tm_rev) = design_cloning_primers(
                        insert_fwd, insert_rev, enz_fwd, enz_rev, primer_conc_nM, custom_prefix=prefix
                    )
                    amplicon_size = len(sequence)

                # Calculate GC content
                gc_fwd = calculate_gc(fwd)
                gc_rev = calculate_gc(rev)

                # Store results
                st.session_state["primer_results"] = {
                    "fwd": fwd, "rev": rev, "len_f": fwd_len_eff, "len_r": rev_len_eff,
                    "gc_f": gc_fwd, "gc_r": gc_rev, "tm_f": tm_fwd, "tm_r": tm_rev,
                    "size": amplicon_size, "mode": mode, "enz_fwd": enz_fwd, 
                    "enz_rev": enz_rev, "prefix": prefix
                }
                
                st.success("✅ Primers successfully generated!")

            except Exception as e:
                st.error(f"❌ Primer design failed: {e}")
                with st.expander("🔍 Debug Information", expanded=False):
                    st.exception(e)
                return
    
    # Display results with all fixes
    if "primer_results" in st.session_state:
        render_primer_results(st.session_state["primer_results"])
    
    # Enhanced usage guide
    st.markdown("---")
    st.markdown("""
    <div class="info-panel">
        <h4>🎯 Design Mode Guide</h4>
        <ul>
            <li><strong>Manual Mode:</strong> Extracts primers directly from sequence termini - fastest option</li>
            <li><strong>Automatic Mode:</strong> AI-optimized primer selection with comprehensive scoring algorithm</li>
            <li><strong>Cloning Mode:</strong> Adds restriction sites and stuffer sequences for molecular cloning</li>
        </ul>
        <p style="margin: 1rem 0 0 0; font-size: 0.9rem; opacity: 0.8;">
            💡 <strong>Tip:</strong> For cloning mode, stuffer sequences are highlighted in <span style="background: #FFD700; padding: 2px 4px; border-radius: 4px;">gold</span> 
            and restriction sites in <span style="background: #FFFF00; padding: 2px 4px; border-radius: 4px;">yellow</span>.
        </p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()