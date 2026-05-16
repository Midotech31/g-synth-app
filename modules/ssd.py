# modules/ssd.py

import re
import logging
from typing import Tuple, Dict, Optional, List

import streamlit as st

from utils.bio_utils import calculate_tm_consensus

logger = logging.getLogger("G‐Synth:SSD")

# ──────────────────────────────────────────────────────────────────────────────
# SSD Constants (exactly as in original)
# ──────────────────────────────────────────────────────────────────────────────

SSD_HIS_TAG = "CACCACCACCACCACCAC"
SSD_LEFT_LINKER = "GGTTCTTCT"
SSD_RIGHT_LINKER = "TCTTCTGGT"

SSD_RESTRICTION_ENZYMES: Dict[str, Dict[str, str]] = {
    "NdeI":    {"recognition": "CATATG",   "cut_forward": "TATG",   "cut_reverse": "CA"},
    "XhoI":    {"recognition": "CTCGAG",   "cut_forward": "C",      "cut_reverse": "TCGAG"},
    "EcoRI":   {"recognition": "GAATTC",   "cut_forward": "AATT",   "cut_reverse": "G"},
    "BamHI":   {"recognition": "GGATCC",   "cut_forward": "GATC",   "cut_reverse": "C"},
    "HindIII": {"recognition": "AAGCTT",   "cut_forward": "AGCT",   "cut_reverse": "T"},
    "SalI":    {"recognition": "GTCGAC",   "cut_forward": "TCGA",   "cut_reverse": "C"},
    "XbaI":    {"recognition": "TCTAGA",   "cut_forward": "CTAG",   "cut_reverse": "A"},
    "NcoI":    {"recognition": "CCATGG",   "cut_forward": "CATG",   "cut_reverse": "G"},
    "KpnI":    {"recognition": "GGTACC",   "cut_forward": "",       "cut_reverse": "GGTAC"},
    "SacI":    {"recognition": "GAGCTC",   "cut_forward": "",       "cut_reverse": "GAGCT"},
    "NotI":    {"recognition": "GCGGCCGC", "cut_forward": "GGCC",   "cut_reverse": "GC"},
    "SpeI":    {"recognition": "ACTAGT",   "cut_forward": "CTAG",   "cut_reverse": "T"},
    "PstI":    {"recognition": "CTGCAG",   "cut_forward": "",       "cut_reverse": "CTGCA"},
    "BglII":   {"recognition": "AGATCT",   "cut_forward": "GATC",   "cut_reverse": "T"},
    "SmaI":    {"recognition": "CCCGGG",   "cut_forward": "",       "cut_reverse": ""},
    "ApaI":    {"recognition": "GGGCCC",   "cut_forward": "",       "cut_reverse": ""},
    "MluI":    {"recognition": "ACGCGT",   "cut_forward": "",       "cut_reverse": ""},
    "EcoRV":   {"recognition": "GATATC",   "cut_forward": "",       "cut_reverse": ""},
    "HpaII":   {"recognition": "CCGG",     "cut_forward": "",       "cut_reverse": ""},
    "SspI":    {"recognition": "AATATT",   "cut_forward": "",       "cut_reverse": ""},
    "DdeI":    {"recognition": "CTNAG",    "cut_forward": "",       "cut_reverse": ""},
    "Bsu36I":  {"recognition": "CCTNAGG",  "cut_forward": "",       "cut_reverse": ""}
}

SSD_CLEAVAGE_SITES: Dict[str, str] = {
    "Thrombin":     "CTGGTGCCGCGTGGTTCT",
    "TEV":          "GAAAACCTGTATTTTCAGGGC",
    "Factor Xa":    "ATCGAAGGTCGT",
    "PreScission":  "CTGGAAGTGCTGTTCCAGGGCCCA",
    "Enterokinase": "GATGACGATGACAAG",
    "SUMO":         "CTGCAGGACTCAGAGG",
    "HRV 3C":       "CTGGAAGTTCTGTTCCAGGGGCCC"
}

# ──────────────────────────────────────────────────────────────────────────────
# Enhanced CSS for SSD Module
# ──────────────────────────────────────────────────────────────────────────────

def inject_ssd_css():
    """Enhanced CSS specifically for SSD module"""
    st.markdown("""
    <style>
    /* SSD-specific enhancements */
    .ssd-input-card {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%) !important;
        border: 2px solid #e2e8f0 !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08) !important;
        transition: all 0.3s ease !important;
    }
    
    .ssd-input-card:hover {
        border-color: #ff4b4b !important;
        box-shadow: 0 8px 25px rgba(255, 75, 75, 0.15) !important;
        transform: translateY(-2px) !important;
    }
    
    .ssd-sequence-display {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%) !important;
        border: 2px solid #ff4b4b !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 8px 25px rgba(255, 75, 75, 0.2) !important;
        font-family: 'Courier New', monospace !important;
        font-size: 14px !important;
        line-height: 1.6 !important;
        color: #ffffff !important;
        overflow-x: auto !important;
        white-space: pre-wrap !important;
        word-break: break-all !important;
    }
    
    .ssd-properties-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.3) !important;
        text-align: center !important;
    }
    
    .ssd-legend-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%) !important;
        border: 2px solid #dee2e6 !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 4px 15px rgba(0,0,0,0.05) !important;
    }
    
    .ssd-section-header {
        background: linear-gradient(135deg, #ff4444 0%, #cc0000 100%) !important;
        color: white !important;
        padding: 0.75rem 1.5rem !important;
        border-radius: 12px !important;
        margin: 1rem 0 0.5rem 0 !important;
        box-shadow: 0 4px 15px rgba(255, 68, 68, 0.3) !important;
        text-align: center !important;
        font-weight: 700 !important;
        font-size: 1.1rem !important;
    }
    
    .ssd-metric-container {
        background: rgba(255, 255, 255, 0.2) !important;
        border-radius: 12px !important;
        padding: 1rem !important;
        margin: 0.5rem 0 !important;
        border: 1px solid rgba(255, 255, 255, 0.3) !important;
        backdrop-filter: blur(10px) !important;
    }
    
    .ssd-highlight-legend {
        display: inline-block;
        padding: 0.3rem 0.8rem;
        border-radius: 8px;
        margin: 0.25rem;
        font-family: 'Courier New', monospace;
        font-weight: 600;
        font-size: 0.9rem;
    }
    
    .ssd-info-box {
        background: linear-gradient(135deg, rgba(34, 197, 94, 0.1), rgba(16, 185, 129, 0.1)) !important;
        border-left: 4px solid #10b981 !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #059669 !important;
    }
    
    .ssd-warning-box {
        background: linear-gradient(135deg, rgba(251, 146, 60, 0.1), rgba(245, 158, 11, 0.1)) !important;
        border-left: 4px solid #f59e0b !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #d97706 !important;
    }
    
    .ssd-options-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 1rem;
        margin: 1rem 0;
    }
    </style>
    """, unsafe_allow_html=True)

# ──────────────────────────────────────────────────────────────────────────────
# Core bio‐sequence helpers (preserved exactly)
# ──────────────────────────────────────────────────────────────────────────────

def ssd_reverse_complement(seq: str) -> str:
    """Reverse complement (handles A/T/C/G/N)."""
    table = str.maketrans("ACGTN", "TGCAN")
    return seq.upper().translate(table)[::-1]

def ssd_validate_sequence(
    sequence: str,
    allow_empty: bool = False,
    allow_ambiguous: bool = False
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate and clean a DNA sequence.
    Returns (is_valid, cleaned_seq, error_message_or_None).
    """
    if not sequence and not allow_empty:
        return False, "", "Sequence cannot be empty"

    valid_chars = "ATCG" + ("RYSWKMBDHVN" if allow_ambiguous else "")
    cleaned = "".join(c for c in sequence.upper() if c in valid_chars)
    if not cleaned and sequence:
        return False, "", "Sequence contains no valid DNA characters"
    if len(cleaned) < len(sequence.replace(" ", "")):
        removed = len(sequence.replace(" ", "")) - len(cleaned)
        return True, cleaned, f"Removed {removed} invalid character{'s' if removed != 1 else ''}"
    return True, cleaned, None

def ssd_calculate_gc_content(seq: str) -> float:
    """Return GC% of a sequence (0.0 if empty)."""
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) * 100 if s else 0.0

# ──────────────────────────────────────────────────────────────────────────────
# Enhanced SSD processing functions
# ──────────────────────────────────────────────────────────────────────────────

def ssd_process_coding_sequence(
    sequence: str,
    remove_stop: bool,
    left_enzyme: str,
    right_enzyme: str
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Process a coding DNA sequence (preserved exactly)"""
    seq = sequence.upper()

    if left_enzyme == "NdeI" and seq.startswith("ATG"):
        seq = seq[3:]

    if remove_stop:
        stop_codons = {"TAA", "TAG", "TGA"}
        modified = seq
        stop_found = False
        for i in range(0, len(modified) - 2, 3):
            codon = modified[i: i + 3]
            if codon in stop_codons:
                modified = modified[:i]
                stop_found = True
                break
        if not stop_found:
            logger.info("No stop codon found to remove.")
        seq = modified

    if left_enzyme not in SSD_RESTRICTION_ENZYMES:
        return None, None, f"Unknown restriction enzyme: {left_enzyme}"
    if right_enzyme not in SSD_RESTRICTION_ENZYMES:
        return None, None, f"Unknown restriction enzyme: {right_enzyme}"

    fwd_cut = SSD_RESTRICTION_ENZYMES[left_enzyme]["cut_forward"]
    rev_cut = SSD_RESTRICTION_ENZYMES[right_enzyme]["cut_forward"]
    fwd_seq = fwd_cut + seq + rev_cut

    rev_cut_r = SSD_RESTRICTION_ENZYMES[right_enzyme]["cut_reverse"]
    fwd_cut_l = SSD_RESTRICTION_ENZYMES[left_enzyme]["cut_reverse"]
    rev_comp = ssd_reverse_complement(seq)
    rev_seq = rev_cut_r + rev_comp + fwd_cut_l

    return fwd_seq, rev_seq, None

def ssd_process_non_coding_sequence(
    sequence: str,
    left_enzyme: str,
    right_enzyme: str,
    cleavage_site: Optional[str],
    include_his_tag: bool = True,
    include_linkers: bool = True
) -> Tuple[str, str]:
    """Enhanced non-coding sequence processing with optional His-tag and linkers"""
    seq = sequence.upper()
    atg_prefix = "" if left_enzyme == "NdeI" else "ATG"
    
    fwd_cut_l = SSD_RESTRICTION_ENZYMES[left_enzyme]["cut_forward"]
    forward = fwd_cut_l + atg_prefix
    
    # Add linkers and His-tag optionally
    if include_linkers and include_his_tag:
        forward += SSD_LEFT_LINKER + SSD_HIS_TAG + SSD_RIGHT_LINKER
    elif include_his_tag:
        # His-tag without linkers
        forward += SSD_HIS_TAG
    elif include_linkers:
        # Linkers without His-tag - use linkers as spacers
        forward += SSD_LEFT_LINKER + SSD_RIGHT_LINKER

    # Add cleavage site if specified
    if cleavage_site and cleavage_site in SSD_CLEAVAGE_SITES:
        forward += SSD_CLEAVAGE_SITES[cleavage_site]

    fwd_cut_r = SSD_RESTRICTION_ENZYMES[right_enzyme]["cut_forward"]
    forward += seq + fwd_cut_r

    # Build reverse strand
    rev_cut_r = SSD_RESTRICTION_ENZYMES[right_enzyme]["cut_reverse"]
    rev_cut_l = SSD_RESTRICTION_ENZYMES[left_enzyme]["cut_reverse"]
    rev = rev_cut_r + ssd_reverse_complement(seq)

    # Add reverse complements in reverse order
    if cleavage_site and cleavage_site in SSD_CLEAVAGE_SITES:
        rev += ssd_reverse_complement(SSD_CLEAVAGE_SITES[cleavage_site])

    if include_linkers and include_his_tag:
        rev += ssd_reverse_complement(SSD_RIGHT_LINKER)
        rev += ssd_reverse_complement(SSD_HIS_TAG)
        rev += ssd_reverse_complement(SSD_LEFT_LINKER)
    elif include_his_tag:
        rev += ssd_reverse_complement(SSD_HIS_TAG)
    elif include_linkers:
        rev += ssd_reverse_complement(SSD_RIGHT_LINKER)
        rev += ssd_reverse_complement(SSD_LEFT_LINKER)

    rev += ssd_reverse_complement(atg_prefix)
    rev += rev_cut_l

    return forward, rev

def ssd_process_sequence(
    input_sequence: str,
    is_coding: bool,
    remove_stop: bool,
    enzyme_pair: str,
    cleavage_site: Optional[str] = None,
    include_his_tag: bool = True,
    include_linkers: bool = True
) -> Dict:
    """Enhanced top‐level entry point for SSD with optional features"""
    valid, seq, err = ssd_validate_sequence(input_sequence)
    if not valid:
        return {"error": err}

    ep = enzyme_pair.replace(" ", "")
    try:
        left_enzyme, right_enzyme = ep.split("/")
    except ValueError:
        return {"error": "Enzyme pair format must be 'Enz1/Enz2'"}

    if left_enzyme not in SSD_RESTRICTION_ENZYMES:
        return {"error": f"Unknown restriction enzyme: {left_enzyme}"}
    if right_enzyme not in SSD_RESTRICTION_ENZYMES:
        return {"error": f"Unknown restriction enzyme: {right_enzyme}"}

    if is_coding:
        fwd, rev, coding_err = ssd_process_coding_sequence(seq, remove_stop, left_enzyme, right_enzyme)
        if coding_err:
            return {"error": coding_err}
    else:
        fwd, rev = ssd_process_non_coding_sequence(
            seq, left_enzyme, right_enzyme, cleavage_site, include_his_tag, include_linkers
        )

    fwd_len = len(fwd)
    rev_len = len(rev)
    fwd_gc = ssd_calculate_gc_content(fwd)
    rev_gc = ssd_calculate_gc_content(rev)

    try:
        fwd_tm = calculate_tm_consensus(fwd)
    except Exception:
        fwd_tm = None

    try:
        rev_tm = calculate_tm_consensus(ssd_reverse_complement(rev))
    except Exception:
        rev_tm = None

    return {
        "forward": fwd,
        "reverse": rev,
        "properties": {
            "forward_length": fwd_len,
            "reverse_length": rev_len,
            "forward_gc": round(fwd_gc, 1),
            "reverse_gc": round(rev_gc, 1),
            "forward_tm": round(fwd_tm, 1) if fwd_tm is not None else None,
            "reverse_tm": round(rev_tm, 1) if rev_tm is not None else None,
        },
    }

# ──────────────────────────────────────────────────────────────────────────────
# Enhanced highlighting functions
# ──────────────────────────────────────────────────────────────────────────────

def highlight_ssd_sequence(seq: str) -> str:
    """Enhanced HTML formatting with improved colors and visibility"""
    highlights: List[Tuple[int, int, str]] = []

    # 1) Restriction sites - Bright yellow with dark text
    for enzyme, props in SSD_RESTRICTION_ENZYMES.items():
        recog = props["recognition"]
        for m in re.finditer(re.escape(recog), seq):
            highlights.append((m.start(), len(recog), "#FFD700", "#000000"))

    # 2) Find first ORF
    orf_start = -1
    orf_stop = -1
    first_atg = seq.find("ATG")
    if first_atg != -1:
        pos = first_atg
        while pos <= len(seq) - 3:
            codon = seq[pos: pos + 3]
            if codon in ("TAA", "TAG", "TGA"):
                orf_start = first_atg
                orf_stop = pos
                break
            pos += 3

        if orf_start != -1:
            highlights.append((orf_start, 3, "#22C55E", "#FFFFFF"))  # Green with white text
        if orf_stop != -1:
            highlights.append((orf_stop, 3, "#EF4444", "#FFFFFF"))  # Red with white text

    # 3) His‐tag - Pink/magenta
    hpos = seq.find(SSD_HIS_TAG)
    if hpos >= 0:
        highlights.append((hpos, len(SSD_HIS_TAG), "#E91E63", "#FFFFFF"))

    # 4) Linkers - Blue
    for linker in (SSD_LEFT_LINKER, SSD_RIGHT_LINKER):
        for m in re.finditer(re.escape(linker), seq):
            highlights.append((m.start(), len(linker), "#2196F3", "#FFFFFF"))

    # 5) Cleavage Sites - Gray
    for site_seq in SSD_CLEAVAGE_SITES.values():
        for m in re.finditer(re.escape(site_seq), seq):
            highlights.append((m.start(), len(site_seq), "#6B7280", "#FFFFFF"))

    # Sort highlights by start position
    highlights.sort(key=lambda x: (x[0], -x[1]))

    html_parts: List[str] = []
    idx = 0
    for start, length, bg_color, text_color in highlights:
        if start < idx:
            continue
        if idx < start:
            html_parts.append(seq[idx:start])
        segment = seq[start: start + length]
        span = f"<span style='background-color:{bg_color}; color:{text_color}; padding:2px 4px; border-radius:4px; font-weight:600;'>{segment}</span>"
        html_parts.append(span)
        idx = start + length

    if idx < len(seq):
        html_parts.append(seq[idx:])

    return "".join(html_parts)

# ──────────────────────────────────────────────────────────────────────────────
# Enhanced Streamlit Interface
# ──────────────────────────────────────────────────────────────────────────────

def render_input_section():
    """Enhanced input section with optional features for non-coding sequences"""
    st.markdown('<div class="ssd-section-header">🔧 Sequence Configuration</div>', unsafe_allow_html=True)
    
    with st.container():
        st.markdown('<div class="ssd-input-card">', unsafe_allow_html=True)
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            seq_type = st.radio(
                "**Sequence Type**",
                ("Coding (starts with ATG)", "Non-Coding"),
                index=0,
                help="Choose whether your sequence codes for a protein or is non-coding"
            )
            
            enzyme_pair = st.selectbox(
                "**Restriction Enzymes**",
                ["NdeI/XhoI", "NdeI/EcoRI", "BamHI/EcoRI", "BamHI/XhoI", "SalI/XbaI"],
                index=0,
                help="Select the restriction enzyme pair for cloning"
            )
        
        with col2:
            is_coding = seq_type.startswith("Coding")
            
            if is_coding:
                remove_stop = st.checkbox(
                    "**Remove in‐frame stop codon**",
                    value=True,
                    help="Automatically remove the first in-frame stop codon"
                )
                cleavage_site = None
                include_his_tag = True
                include_linkers = True
                
                if remove_stop:
                    st.markdown('<div class="ssd-info-box">ℹ️ First in-frame stop codon will be removed automatically</div>', unsafe_allow_html=True)
            else:
                remove_stop = False
                
                # Non-coding sequence options
                st.markdown("**Expression Features:**")
                
                # His-tag option
                include_his_tag = st.checkbox(
                    "Include His-Tag",
                    value=True,
                    help="Add 6xHistidine tag for protein purification"
                )
                
                # Linkers option
                include_linkers = st.checkbox(
                    "Include Flexible Linkers",
                    value=True,
                    help="Add flexible linker sequences around His-tag"
                )
                
                # Cleavage site option
                cleavage_options = ["None"] + list(SSD_CLEAVAGE_SITES.keys())
                cleavage_selection = st.selectbox(
                    "**Protease Cleavage Site**",
                    cleavage_options,
                    index=0,
                    help="Choose a protease cleavage site for protein purification, or None to skip"
                )
                cleavage_site = None if cleavage_selection == "None" else cleavage_selection
                
                # Show status messages
                if not include_his_tag and not include_linkers and cleavage_site is None:
                    st.markdown('<div class="ssd-warning-box">⚠️ No expression features selected - sequence will only have restriction sites</div>', unsafe_allow_html=True)
                else:
                    features = []
                    if include_his_tag:
                        features.append("His-Tag")
                    if include_linkers:
                        features.append("Linkers")
                    if cleavage_site:
                        features.append(f"{cleavage_site} cleavage site")
                    
                    if features:
                        features_text = ", ".join(features)
                        st.markdown(f'<div class="ssd-info-box">✅ Including: {features_text}</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Sequence input area
        st.markdown('<div class="ssd-input-card">', unsafe_allow_html=True)
        seq_input = st.text_area(
            "**Enter DNA Sequence (5' → 3')**",
            value="",
            height=120,
            help="Enter your DNA sequence using only A, T, C, G characters",
            placeholder="ATGAAACCGGATTTACGTCAA..."
        )
        st.markdown('</div>', unsafe_allow_html=True)
        
        return seq_input, is_coding, remove_stop, enzyme_pair, cleavage_site, include_his_tag, include_linkers

def render_sequence_results(fwd_seq: str, rev_seq: str):
    """Enhanced sequence results display"""
    st.markdown('<div class="ssd-section-header">🧬 Generated Sequences</div>', unsafe_allow_html=True)
    
    # Forward strand
    st.markdown("**Forward Strand (5′ → 3′)**")
    highlighted_fwd = highlight_ssd_sequence(fwd_seq)
    st.markdown(f'<div class="ssd-sequence-display">{highlighted_fwd}</div>', unsafe_allow_html=True)
    
    # Add copy button for forward sequence
    if st.button("📋 Copy Forward Sequence", key="copy_fwd"):
        st.code(fwd_seq, language="text")
    
    # Reverse strand
    st.markdown("**Reverse Strand (5′ → 3′)**")
    highlighted_rev = highlight_ssd_sequence(rev_seq)
    st.markdown(f'<div class="ssd-sequence-display">{highlighted_rev}</div>', unsafe_allow_html=True)
    
    # Add copy button for reverse sequence
    if st.button("📋 Copy Reverse Sequence", key="copy_rev"):
        st.code(rev_seq, language="text")

def render_properties_section(props: Dict):
    """Enhanced properties display with modern styling"""
    st.markdown('<div class="ssd-section-header">📊 Sequence Properties</div>', unsafe_allow_html=True)
    
    st.markdown('<div class="ssd-properties-card">', unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown('<div class="ssd-metric-container">', unsafe_allow_html=True)
        st.markdown("**Forward Strand**")
        st.markdown(f"**Length:** {props['forward_length']} bp")
        st.markdown(f"**GC Content:** {props['forward_gc']:.1f}%")
        fwtm = f"{props['forward_tm']:.1f}°C" if props["forward_tm"] is not None else "N/A"
        st.markdown(f"**Melting Temp:** {fwtm}")
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col2:
        st.markdown('<div class="ssd-metric-container">', unsafe_allow_html=True)
        st.markdown("**Reverse Strand**")
        st.markdown(f"**Length:** {props['reverse_length']} bp")
        st.markdown(f"**GC Content:** {props['reverse_gc']:.1f}%")
        rvtm = f"{props['reverse_tm']:.1f}°C" if props["reverse_tm"] is not None else "N/A"
        st.markdown(f"**Melting Temp:** {rvtm}")
        st.markdown('</div>', unsafe_allow_html=True)
        
    with col3:
        st.markdown('<div class="ssd-metric-container">', unsafe_allow_html=True)
        st.markdown("**Quality Metrics**")
        
        # GC content assessment
        avg_gc = (props['forward_gc'] + props['reverse_gc']) / 2
        if 40 <= avg_gc <= 60:
            gc_status = "🟢 Optimal"
        elif 30 <= avg_gc <= 70:
            gc_status = "🟡 Acceptable"
        else:
            gc_status = "🔴 Suboptimal"
        
        st.markdown(f"**GC Status:** {gc_status}")
        st.markdown(f"**Average GC:** {avg_gc:.1f}%")
        
        # Length assessment
        if props['forward_length'] < 100:
            length_status = "🟢 Short"
        elif props['forward_length'] < 500:
            length_status = "🟡 Medium"
        else:
            length_status = "🔴 Long"
        
        st.markdown(f"**Length:** {length_status}")
        st.markdown('</div>', unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)

def render_enhanced_legend():
    """Enhanced color legend with modern styling"""
    st.markdown('<div class="ssd-section-header">🎨 Sequence Highlighting Legend</div>', unsafe_allow_html=True)
    
    st.markdown('<div class="ssd-legend-card">', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Functional Elements:**")
        st.markdown('<span class="ssd-highlight-legend" style="background-color:#22C55E; color:#FFFFFF;">ATG</span> Start codon (first ORF)', unsafe_allow_html=True)
        st.markdown('<span class="ssd-highlight-legend" style="background-color:#EF4444; color:#FFFFFF;">TAA</span> Stop codon (first ORF)', unsafe_allow_html=True)
        st.markdown('<span class="ssd-highlight-legend" style="background-color:#E91E63; color:#FFFFFF;">CACCAC...</span> His-Tag', unsafe_allow_html=True)
        
    with col2:
        st.markdown("**Cloning Elements:**")
        st.markdown('<span class="ssd-highlight-legend" style="background-color:#FFD700; color:#000000;">GAATTC</span> Restriction sites', unsafe_allow_html=True)
        st.markdown('<span class="ssd-highlight-legend" style="background-color:#2196F3; color:#FFFFFF;">GGTTCTTCT</span> Linkers', unsafe_allow_html=True)
        st.markdown('<span class="ssd-highlight-legend" style="background-color:#6B7280; color:#FFFFFF;">CTGGAA...</span> Cleavage sites', unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)

def main():
    """Enhanced Streamlit entry‐point with modern design"""
    # Inject enhanced CSS
    inject_ssd_css()
    
    # Main header with description
    st.markdown("""
    <div style="
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 2rem;
        border-radius: 16px;
        margin-bottom: 2rem;
        text-align: center;
        box-shadow: 0 8px 32px rgba(102, 126, 234, 0.3);
    ">
        <h2 style="margin: 0; font-size: 2rem; font-weight: 800;">
            🔬 Small Sequence Design
        </h2>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.1rem; opacity: 0.9;">
            Design optimized DNA sequences for molecular cloning with customizable expression features
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Input section
    seq_input, is_coding, remove_stop, enzyme_pair, cleavage_site, include_his_tag, include_linkers = render_input_section()
    
    # Process button with enhanced styling
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        process_clicked = st.button(
            "🚀 Process Sequence",
            use_container_width=True,
            type="primary"
        )
    
    if process_clicked:
        raw = seq_input.strip().upper()
        if not raw:
            st.error("Please enter a DNA sequence.")
            st.stop()

        clean_seq = re.sub(r"[^ATCG]", "", raw)
        if not clean_seq:
            st.error("Invalid DNA sequence. Use only A, T, C, G.")
            st.stop()

        if is_coding and not clean_seq.startswith("ATG"):
            st.error("Coding sequence must start with ATG.")
            st.stop()

        # Show processing indicator
        with st.spinner("Processing sequence..."):
            result = ssd_process_sequence(
                input_sequence=clean_seq,
                is_coding=is_coding,
                remove_stop=remove_stop,
                enzyme_pair=enzyme_pair,
                cleavage_site=cleavage_site,
                include_his_tag=include_his_tag,
                include_linkers=include_linkers
            )

        if "error" in result:
            st.error(result["error"])
            st.stop()

        # Display results with enhanced styling
        fwd_seq = result["forward"]
        rev_seq = result["reverse"]
        props = result["properties"]
        
        render_sequence_results(fwd_seq, rev_seq)
        render_properties_section(props)
        
        # Add download functionality
        st.markdown('<div class="ssd-section-header">💾 Export Options</div>', unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        with col1:
            # Format sequences for download
            fasta_content = f">Forward_Strand_5to3\n{fwd_seq}\n>Reverse_Strand_5to3\n{rev_seq}"
            st.download_button(
                "📄 Download FASTA",
                fasta_content,
                file_name="ssd_sequences.fasta",
                mime="text/plain"
            )
        
        with col2:
            # Create enhanced summary report
            features_used = []
            if include_his_tag:
                features_used.append("His-Tag")
            if include_linkers:
                features_used.append("Flexible Linkers")
            if cleavage_site:
                features_used.append(f"{cleavage_site} cleavage site")
            
            features_text = ", ".join(features_used) if features_used else "None"
            
            summary = f"""Small Sequence Design Report
{'='*50}
Input Sequence: {clean_seq}
Sequence Type: {'Coding' if is_coding else 'Non-coding'}
Enzyme Pair: {enzyme_pair}
Features Used: {features_text}

Forward Strand (5' -> 3'): {fwd_seq}
Reverse Strand (5' -> 3'): {rev_seq}

Properties:
- Forward: {props['forward_length']} bp, {props['forward_gc']:.1f}% GC
- Reverse: {props['reverse_length']} bp, {props['reverse_gc']:.1f}% GC
"""
            st.download_button(
                "📊 Download Report",
                summary,
                file_name="ssd_report.txt",
                mime="text/plain"
            )
    
    # Enhanced legend
    render_enhanced_legend()
    
    # Additional help section
    with st.expander("💡 Usage Tips & Best Practices", expanded=False):
        st.markdown("""
        **For Coding Sequences:**
        - Must start with ATG (start codon)
        - Automatically handles stop codon removal
        - NdeI enzyme pair automatically removes initial ATG
        
        **For Non-Coding Sequences:**
        - **His-Tag**: Adds 6xHistidine tag for protein purification
        - **Flexible Linkers**: Add spacing around His-tag for better accessibility
        - **Cleavage Sites**: Enable tag removal with specific proteases
        - **Minimal Design**: Uncheck all features for basic cloning with restriction sites only
        
        **Quality Guidelines:**
        - Optimal GC content: 40-60%
        - Avoid long stretches of identical nucleotides
        - Check for unwanted restriction sites in your sequence
        """)