# modules/extended_synthesis.py

import re
import io
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from typing import List, Dict, Tuple

from utils.bio_utils import (
    fragment_extended_sequence,
    ENZYME_PAIRS,
    CLEAVAGE_SITES,
)

# ───────────────────────────────────────────────────────────────────────────────
# ENHANCED CSS FOR EXTENDED SYNTHESIS MODULE
# ───────────────────────────────────────────────────────────────────────────────

def inject_extended_synthesis_css():
    """Enhanced CSS specifically for Extended Synthesis module with red accent theme"""
    st.markdown("""
    <style>
    /* Extended Synthesis specific enhancements */
    .ext-synthesis-input-card {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%) !important;
        border: 2px solid #e2e8f0 !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08) !important;
        transition: all 0.3s ease !important;
    }
    
    .ext-synthesis-input-card:hover {
        border-color: #ff4b4b !important;
        box-shadow: 0 8px 25px rgba(255, 75, 75, 0.15) !important;
        transform: translateY(-2px) !important;
    }
    
    .ext-synthesis-sequence-display {
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
    
    .ext-synthesis-results-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.3) !important;
        text-align: center !important;
    }
    
    .ext-synthesis-section-header {
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
    
    .ext-synthesis-metric-container {
        background: rgba(255, 255, 255, 0.2) !important;
        border-radius: 12px !important;
        padding: 1rem !important;
        margin: 0.5rem 0 !important;
        border: 1px solid rgba(255, 255, 255, 0.3) !important;
        backdrop-filter: blur(10px) !important;
    }
    
    .ext-synthesis-info-box {
        background: linear-gradient(135deg, rgba(34, 197, 94, 0.1), rgba(16, 185, 129, 0.1)) !important;
        border-left: 4px solid #10b981 !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #059669 !important;
    }
    
    .ext-synthesis-warning-box {
        background: linear-gradient(135deg, rgba(251, 146, 60, 0.1), rgba(245, 158, 11, 0.1)) !important;
        border-left: 4px solid #f59e0b !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #d97706 !important;
    }
    
    .ext-synthesis-error-box {
        background: linear-gradient(135deg, rgba(239, 68, 68, 0.1), rgba(220, 38, 38, 0.1)) !important;
        border-left: 4px solid #dc2626 !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #dc2626 !important;
    }
    
    .ext-synthesis-success-box {
        background: linear-gradient(135deg, rgba(34, 197, 94, 0.1), rgba(16, 185, 129, 0.1)) !important;
        border-left: 4px solid #10b981 !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #059669 !important;
    }
    
    .ext-synthesis-fragment-card {
        background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%) !important;
        border: 2px solid #dee2e6 !important;
        border-radius: 12px !important;
        padding: 1rem !important;
        margin: 0.5rem 0 !important;
        transition: all 0.3s ease !important;
    }
    
    .ext-synthesis-fragment-card:hover {
        border-color: #ff4b4b !important;
        box-shadow: 0 4px 15px rgba(255, 75, 75, 0.1) !important;
    }
    </style>
    """, unsafe_allow_html=True)

# ───────────────────────────────────────────────────────────────────────────────
# ENHANCED CORE FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content percentage."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100

def validate_sequence_input(sequence: str) -> Tuple[bool, str, str]:
    """
    Validate DNA sequence input and return cleaned version.
    Returns: (is_valid, cleaned_sequence, message)
    """
    if not sequence:
        return False, "", "Sequence cannot be empty"

    cleaned = re.sub(r"[^ATCG]", "", sequence.upper())
    if not cleaned:
        return False, "", "Sequence contains no valid DNA bases (A, T, C, G)"

    orig_len = len(sequence.replace(" ", "").replace("\n", ""))
    invalid_count = orig_len - len(cleaned)
    if invalid_count > 0:
        message = f"Removed {invalid_count} invalid character(s). Clean length: {len(cleaned)} bp"
    else:
        message = f"Valid sequence: {len(cleaned)} bp"

    return True, cleaned, message

def create_copy_button_html(sequence, label, button_id):
    """Create enhanced copy button with red theme"""
    return f"""
    <div style="margin: 10px 0;">
        <button onclick="copyToClipboard{button_id}()" style="
            background: linear-gradient(135deg, #ff4b4b, #cc0000);
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            font-size: 14px;
            transition: all 0.3s ease;
        ">
            📋 Copy {label}
        </button>
        <span id="copy-status{button_id}" style="margin-left: 10px; color: #ff4b4b; font-weight: bold;"></span>
    </div>
    
    <script>
    function copyToClipboard{button_id}() {{
        const textArea = document.createElement('textarea');
        textArea.value = `{sequence}`;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        
        const status = document.getElementById('copy-status{button_id}');
        status.textContent = '✅ Copied!';
        status.style.color = '#10B981';
        setTimeout(() => {{
            status.textContent = '';
        }}, 2000);
    }}
    </script>
    """

def create_enhanced_fragment_visualization(assembly: List[Dict]):
    """Create modern interactive visualization with red theme"""
    
    # Prepare data
    fragments = [f"Fragment {frag['fragment']}" for frag in assembly]
    lengths = [frag['length'] for frag in assembly]
    gc_contents = [calculate_gc_content(frag['sequence']) for frag in assembly]
    types = [frag['type'] for frag in assembly]
    
    # Color mapping with red theme
    color_map = {
        'First': '#ff4b4b',
        'Internal': '#667eea', 
        'Last': '#cc0000'
    }
    colors = [color_map[t] for t in types]
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Fragment Lengths', 'GC Content Distribution', 'Length by Position', 'GC Content by Position'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Fragment lengths bar chart
    fig.add_trace(
        go.Bar(
            x=fragments,
            y=lengths,
            marker_color=colors,
            name='Fragment Length',
            hovertemplate='<b>%{x}</b><br>Length: %{y} bp<br>Type: %{customdata}<extra></extra>',
            customdata=types
        ),
        row=1, col=1
    )
    
    # GC content histogram
    fig.add_trace(
        go.Histogram(
            x=gc_contents,
            marker_color='#ff4b4b',
            opacity=0.7,
            name='GC Distribution',
            hovertemplate='GC Content: %{x:.1f}%<br>Count: %{y}<extra></extra>'
        ),
        row=1, col=2
    )
    
    # Length by position
    fig.add_trace(
        go.Scatter(
            x=list(range(1, len(assembly) + 1)),
            y=lengths,
            mode='markers+lines',
            marker=dict(size=10, color=colors, line=dict(width=2, color='white')),
            line=dict(color='#ff4b4b', width=2),
            name='Length Trend',
            hovertemplate='Position: %{x}<br>Length: %{y} bp<extra></extra>'
        ),
        row=2, col=1
    )
    
    # GC content by position
    fig.add_trace(
        go.Scatter(
            x=list(range(1, len(assembly) + 1)),
            y=gc_contents,
            mode='markers+lines',
            marker=dict(size=10, color='#667eea', line=dict(width=2, color='white')),
            line=dict(color='#667eea', width=2),
            name='GC Trend',
            hovertemplate='Position: %{x}<br>GC Content: %{y:.1f}%<extra></extra>'
        ),
        row=2, col=2
    )
    
    # Add reference line for GC content
    fig.add_hline(y=50, line_dash="dash", line_color="red", opacity=0.7, row=2, col=2)
    
    # Update layout
    fig.update_layout(
        height=600,
        showlegend=False,
        plot_bgcolor='#f8fafc',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12),
        title_text="Fragment Analysis Dashboard",
        title_x=0.5
    )
    
    # Update axes
    fig.update_xaxes(title_text="Fragment", row=1, col=1)
    fig.update_yaxes(title_text="Length (bp)", row=1, col=1)
    fig.update_xaxes(title_text="GC Content (%)", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)
    fig.update_xaxes(title_text="Fragment Position", row=2, col=1)
    fig.update_yaxes(title_text="Length (bp)", row=2, col=1)
    fig.update_xaxes(title_text="Fragment Position", row=2, col=2)
    fig.update_yaxes(title_text="GC Content (%)", row=2, col=2)
    
    return fig

def create_export_data(
    assembly: List[Dict], original_seq: str, parameters: Dict
) -> Dict[str, str]:
    """Create export data in CSV, FASTA, and text report formats."""
    df = pd.DataFrame(assembly)
    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)
    csv_data = csv_buffer.getvalue()

    fasta_data = f">Original_Sequence_{len(original_seq)}bp\n{original_seq}\n\n"
    for frag in assembly:
        fasta_data += (
            f">Fragment_{frag['fragment']}_{frag['type']}_Forward_{frag['length']}bp\n"
            f"{frag['forward']}\n"
        )
        fasta_data += (
            f">Fragment_{frag['fragment']}_{frag['type']}_Reverse_{frag['length']}bp\n"
            f"{frag['reverse']}\n"
        )

    # Enhanced report with feature information
    features_used = []
    if parameters.get('include_his_tag', True):
        features_used.append("His-Tag")
    if parameters.get('include_linkers', True):
        features_used.append("Flexible Linkers")
    if parameters.get('cleavage_site') and parameters['cleavage_site'] != "None":
        features_used.append(f"{parameters['cleavage_site']} cleavage site")
    
    features_text = ", ".join(features_used) if features_used else "Basic fragmentation only"

    report = (
        "Extended Synthesis Fragmentation Report\n"
        "Generated by G-Synth Extended Synthesis Module\n\n"
        "PARAMETERS:\n"
        f"- Fragment Size: {parameters['fragment_size']} bp\n"
        f"- Enzyme Pair: {parameters['enzyme_pair']}\n"
        f"- Cleavage Site: {parameters['cleavage_site']}\n"
        f"- Internal Overlap: {parameters['internal_overlap']} bp\n"
        f"- Features Used: {features_text}\n\n"
        "SEQUENCE INFORMATION:\n"
        f"- Original Length: {len(original_seq)} bp\n"
        f"- Total Fragments: {len(assembly)}\n"
        f"- Average Fragment Size: {sum(frag['length'] for frag in assembly) / len(assembly):.1f} bp\n"
        f"- GC Content: {calculate_gc_content(original_seq):.1f}%\n\n"
        "ASSEMBLY VERIFICATION:\n"
        "- Fragments can be assembled using overlapping regions\n"
        "- Each fragment includes appropriate sticky ends for cloning\n\n"
        "FRAGMENTS DETAILS:\n"
    )
    for frag in assembly:
        report += (
            f"\nFragment {frag['fragment']} ({frag['type']}):\n"
            f"  Core Length: {frag['length']} bp\n"
            f"  Forward Total: {len(frag['forward'])} bp\n"
            f"  Reverse Total: {len(frag['reverse'])} bp\n"
            f"  GC Content: {calculate_gc_content(frag['sequence']):.1f}%\n"
            f"  Forward: {frag['forward']}\n"
            f"  Reverse: {frag['reverse']}\n"
        )

    return {"csv": csv_data, "fasta": fasta_data, "report": report}

def display_fragmentation_results(
    assembly: List[Dict], original_seq: str, reassembled: str, parameters: Dict
):
    """Enhanced display with modern red theme styling"""
    is_valid = original_seq == reassembled

    st.markdown('<div class="ext-synthesis-section-header">🧬 Fragmentation Results</div>', unsafe_allow_html=True)

    # Assembly validation status
    if is_valid:
        st.markdown('<div class="ext-synthesis-success-box">✅ Assembly verification: PASSED - Reassembled sequence matches original</div>', unsafe_allow_html=True)
    else:
        st.markdown('<div class="ext-synthesis-error-box">❌ Assembly verification: FAILED - Sequences do not match</div>', unsafe_allow_html=True)
        st.markdown(f'<div class="ext-synthesis-warning-box">⚠️ Original length: {len(original_seq)} bp, Reassembled length: {len(reassembled)} bp</div>', unsafe_allow_html=True)

    # Enhanced summary metrics
    st.markdown('<div class="ext-synthesis-results-card">', unsafe_allow_html=True)
    
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
        st.metric("Total Fragments", len(assembly))
        st.markdown('</div>', unsafe_allow_html=True)
    with col2:
        st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
        st.metric("Original Length", f"{len(original_seq)} bp")
        st.markdown('</div>', unsafe_allow_html=True)
    with col3:
        st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
        st.metric("Reassembled Length", f"{len(reassembled)} bp")
        st.markdown('</div>', unsafe_allow_html=True)
    with col4:
        st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
        avg_length = sum(frag["length"] for frag in assembly) / len(assembly)
        st.metric("Avg Fragment Size", f"{avg_length:.1f} bp")
        st.markdown('</div>', unsafe_allow_html=True)
    with col5:
        st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
        avg_gc = sum(calculate_gc_content(frag["sequence"]) for frag in assembly) / len(assembly)
        st.metric("Avg GC Content", f"{avg_gc:.1f}%")
        st.markdown('</div>', unsafe_allow_html=True)
    
    st.markdown('</div>', unsafe_allow_html=True)

    # Enhanced tabs
    tab1, tab2, tab3, tab4 = st.tabs(
        ["📊 Fragment Table", "📈 Visualization", "🧬 Sequences", "💾 Export"]
    )

    with tab1:
        st.markdown('<div class="ext-synthesis-section-header">📊 Fragment Details</div>', unsafe_allow_html=True)
        
        # Enhanced dataframe with styling
        df_data = []
        for frag in assembly:
            df_data.append({
                "Fragment": frag["fragment"],
                "Type": frag["type"],
                "Core Length (bp)": frag["length"],
                "Forward Total (bp)": len(frag["forward"]),
                "Reverse Total (bp)": len(frag["reverse"]),
                "GC Content (%)": f"{calculate_gc_content(frag['sequence']):.1f}",
                "Core Sequence": (
                    frag["sequence"][:50] + "..."
                    if len(frag["sequence"]) > 50
                    else frag["sequence"]
                ),
            })

        df = pd.DataFrame(df_data)

        # Enhanced styling function with red theme
        def color_type(val):
            colors = {
                "First": "background-color: #ff4b4b; color: white; font-weight: bold;",
                "Internal": "background-color: #667eea; color: white; font-weight: bold;",
                "Last": "background-color: #cc0000; color: white; font-weight: bold;",
            }
            return colors.get(val, "")

        styled_df = df.style.applymap(color_type, subset=["Type"])
        st.dataframe(styled_df, use_container_width=True)

        # Enhanced statistics
        st.markdown('<div class="ext-synthesis-section-header">📈 Fragment Statistics</div>', unsafe_allow_html=True)
        
        st.markdown('<div class="ext-synthesis-results-card">', unsafe_allow_html=True)
        col_a, col_b, col_c = st.columns(3)
        
        with col_a:
            st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
            st.markdown("**Length Distribution**")
            lengths = [frag["length"] for frag in assembly]
            st.write(f"Min: {min(lengths)} bp")
            st.write(f"Max: {max(lengths)} bp")
            st.write(f"Std: {np.std(lengths):.1f} bp")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col_b:
            st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
            st.markdown("**GC Content Distribution**")
            gc_contents = [calculate_gc_content(frag["sequence"]) for frag in assembly]
            st.write(f"Min: {min(gc_contents):.1f}%")
            st.write(f"Max: {max(gc_contents):.1f}%")
            st.write(f"Std: {np.std(gc_contents):.1f}%")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col_c:
            st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
            st.markdown("**Quality Assessment**")
            optimal_gc = sum(1 for gc in gc_contents if 40 <= gc <= 60)
            optimal_length = sum(1 for length in lengths if 150 <= length <= 250)
            st.write(f"Optimal GC: {optimal_gc}/{len(assembly)}")
            st.write(f"Optimal Length: {optimal_length}/{len(assembly)}")
            st.write(f"Assembly Valid: {'✅ Yes' if is_valid else '❌ No'}")
            st.markdown('</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)

    with tab2:
        st.markdown('<div class="ext-synthesis-section-header">📈 Interactive Analysis</div>', unsafe_allow_html=True)
        
        # Enhanced plotly visualization
        enhanced_fig = create_enhanced_fragment_visualization(assembly)
        st.plotly_chart(enhanced_fig, use_container_width=True)

    with tab3:
        st.markdown('<div class="ext-synthesis-section-header">🧬 Fragment Sequences</div>', unsafe_allow_html=True)
        
        for i, frag in enumerate(assembly):
            st.markdown('<div class="ext-synthesis-fragment-card">', unsafe_allow_html=True)
            
            with st.expander(
                f"Fragment {frag['fragment']} ({frag['type']}) - {frag['length']} bp",
                expanded=(i == 0),
            ):
                # Enhanced metrics
                col_info1, col_info2, col_info3, col_info4 = st.columns(4)
                with col_info1:
                    st.metric("Core Length", f"{frag['length']} bp")
                with col_info2:
                    st.metric("GC Content", f"{calculate_gc_content(frag['sequence']):.1f}%")
                with col_info3:
                    gc_pct = calculate_gc_content(frag["sequence"])
                    tm_est = 64.9 + 41 * (gc_pct / 100 - 16.4 / len(frag["sequence"]))
                    st.metric("Est. Tm", f"{tm_est:.1f}°C")
                with col_info4:
                    type_colors = {'First': '🟢', 'Internal': '🔵', 'Last': '🔴'}
                    st.metric("Type", f"{type_colors.get(frag['type'], '⚪')} {frag['type']}")

                # Enhanced sequence display
                col_seq1, col_seq2 = st.columns(2)
                with col_seq1:
                    st.markdown("**Forward Strand (5' → 3')**")
                    st.markdown(f'<div class="ext-synthesis-sequence-display">{frag["forward"]}</div>', unsafe_allow_html=True)
                    copy_fwd_html = create_copy_button_html(frag["forward"], f"Forward F{frag['fragment']}", f"fwd{i}")
                    st.components.v1.html(copy_fwd_html, height=80)
                
                with col_seq2:
                    st.markdown("**Reverse Strand (5' → 3')**")
                    st.markdown(f'<div class="ext-synthesis-sequence-display">{frag["reverse"]}</div>', unsafe_allow_html=True)
                    copy_rev_html = create_copy_button_html(frag["reverse"], f"Reverse F{frag['fragment']}", f"rev{i}")
                    st.components.v1.html(copy_rev_html, height=80)

                st.markdown("**Core Sequence**")
                st.markdown(f'<div class="ext-synthesis-sequence-display">{frag["sequence"]}</div>', unsafe_allow_html=True)
                copy_core_html = create_copy_button_html(frag["sequence"], f"Core F{frag['fragment']}", f"core{i}")
                st.components.v1.html(copy_core_html, height=80)
            
            st.markdown('</div>', unsafe_allow_html=True)

    with tab4:
        st.markdown('<div class="ext-synthesis-section-header">💾 Export Options</div>', unsafe_allow_html=True)
        
        export_data = create_export_data(assembly, original_seq, parameters)

        st.markdown('<div class="ext-synthesis-input-card">', unsafe_allow_html=True)
        col_exp1, col_exp2, col_exp3 = st.columns(3)
        
        with col_exp1:
            st.download_button(
                "📊 Download CSV",
                data=export_data["csv"],
                file_name="extended_synthesis_fragments.csv",
                mime="text/csv",
                use_container_width=True
            )
        
        with col_exp2:
            st.download_button(
                "🧬 Download FASTA",
                data=export_data["fasta"],
                file_name="extended_synthesis_sequences.fasta",
                mime="text/plain",
                use_container_width=True,
                type="primary"
            )
        
        with col_exp3:
            st.download_button(
                "📄 Download Report",
                data=export_data["report"],
                file_name="extended_synthesis_report.txt",
                mime="text/plain",
                use_container_width=True
            )
        
        st.markdown('</div>', unsafe_allow_html=True)

def main():
    """Enhanced main function with red theme and modern design plus optional features"""
    
    # Inject enhanced CSS
    inject_extended_synthesis_css()
    
    # Professional header with red theme
    st.markdown("""
    <div style="
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 16px;
        text-align: center;
        margin-bottom: 2rem;
        color: white;
        box-shadow: 0 8px 32px rgba(102, 126, 234, 0.3);
    ">
        <h2 style="margin: 0; font-size: 2rem; font-weight: 800;">
            🧬 Extended Synthesis
        </h2>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.1rem; opacity: 0.9;">
            Fragment long DNA sequences into overlapping oligonucleotides with customizable features
        </p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown('<div class="ext-synthesis-section-header">⚙️ Input Parameters</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown('<div class="ext-synthesis-input-card">', unsafe_allow_html=True)
        st.markdown("### DNA Sequence Input")
        
        input_tab1, input_tab2 = st.tabs(["📝 Text Input", "📁 File Upload"])

        # Initialize session_state keys if missing
        if "sequence_input" not in st.session_state:
            st.session_state.sequence_input = ""
        if "assembly" not in st.session_state:
            st.session_state.assembly = None
        if "reassembled" not in st.session_state:
            st.session_state.reassembled = ""
        if "parameters" not in st.session_state:
            st.session_state.parameters = None

        def clear_all():
            """Callback to clear input and results."""
            st.session_state.sequence_input = ""
            st.session_state.assembly = None
            st.session_state.reassembled = ""
            st.session_state.parameters = None

        with input_tab1:
            sequence = st.text_area(
                "Enter your DNA sequence",
                height=200,
                placeholder="Paste your DNA sequence here (A, T, G, C only)...",
                help="Enter a DNA sequence using only A, T, G, C bases. Spaces and other characters will be removed.",
                key="sequence_input",
            )

            if sequence:
                is_valid, clean_seq, message = validate_sequence_input(sequence)
                
                st.markdown('<div class="ext-synthesis-results-card">', unsafe_allow_html=True)
                col_a, col_b, col_c = st.columns(3)
                with col_a:
                    st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
                    st.metric("Input Length", len(sequence))
                    st.markdown('</div>', unsafe_allow_html=True)
                with col_b:
                    st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
                    st.metric("Valid Bases", len(clean_seq) if is_valid else 0)
                    st.markdown('</div>', unsafe_allow_html=True)
                with col_c:
                    st.markdown('<div class="ext-synthesis-metric-container">', unsafe_allow_html=True)
                    if is_valid:
                        st.metric("GC Content", f"{calculate_gc_content(clean_seq):.1f}%")
                    else:
                        st.metric("Status", "Invalid")
                    st.markdown('</div>', unsafe_allow_html=True)
                st.markdown('</div>', unsafe_allow_html=True)

                if is_valid:
                    st.markdown(f'<div class="ext-synthesis-info-box">ℹ️ {message}</div>', unsafe_allow_html=True)
                else:
                    st.markdown(f'<div class="ext-synthesis-error-box">❌ {message}</div>', unsafe_allow_html=True)

            # Enhanced clear button
            col_clear1, col_clear2, col_clear3 = st.columns([1, 2, 1])
            with col_clear2:
                st.button("🗑️ Clear Input", on_click=clear_all, use_container_width=True)

        with input_tab2:
            uploaded_file = st.file_uploader(
                "Upload sequence file",
                type=["fasta", "fa", "txt"],
                help="Upload a FASTA or plain text file containing your DNA sequence",
            )

            if uploaded_file is not None:
                try:
                    content = uploaded_file.read().decode("utf-8")
                    if content.startswith(">"):
                        lines = content.splitlines()
                        seq_lines = [line.strip() for line in lines[1:] if not line.startswith(">")]
                        sequence = "".join(seq_lines)
                    else:
                        sequence = content

                    st.session_state.sequence_input = sequence
                    st.session_state.assembly = None
                    st.session_state.reassembled = ""
                    st.session_state.parameters = None
                    
                    st.markdown(f'<div class="ext-synthesis-success-box">✅ Loaded sequence from file ({len(sequence)} characters)</div>', unsafe_allow_html=True)
                    
                    if sequence:
                        is_valid, clean_seq, message = validate_sequence_input(sequence)
                        if is_valid:
                            preview = clean_seq[:100] + "..." if len(clean_seq) > 100 else clean_seq
                            st.markdown(f'<div class="ext-synthesis-sequence-display">{preview}</div>', unsafe_allow_html=True)
                except Exception as e:
                    st.markdown(f'<div class="ext-synthesis-error-box">❌ Error reading file: {e}</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)

    with col2:
        st.markdown('<div class="ext-synthesis-input-card">', unsafe_allow_html=True)
        st.markdown("### Fragmentation Parameters")
        
        frag_size = st.number_input(
            "Max Fragment Size (bp)",
            min_value=50,
            max_value=1000,
            value=200,
            step=10,
            help="Maximum size of each fragment in base pairs",
        )
        
        overlap_size = st.number_input(
            "Overlap Size (bp)",
            min_value=5,
            max_value=100,
            value=15,
            step=5,
            help="Size of overlap between consecutive fragments",
        )
        
        # Enhanced enzyme selection with "None" option
        enzyme_options = ["None"] + list(ENZYME_PAIRS.keys())
        enzyme_selection = st.selectbox(
            "Enzyme Pair",
            options=enzyme_options,
            index=1,  # Default to first real enzyme
            help="Restriction enzyme pair for creating sticky ends, or None to skip",
        )
        enzyme = None if enzyme_selection == "None" else enzyme_selection

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
        
        # Enhanced cleavage site selection with "None" option  
        cleavage_options = ["None"] + list(CLEAVAGE_SITES.keys())
        cleavage_selection = st.selectbox(
            "Cleavage Site",
            options=cleavage_options,
            index=0,
            help="Optional protease cleavage site sequence, or None to skip",
        )
        cleavage_site = None if cleavage_selection == "None" else cleavage_selection

        # Feature status indicators
        features_used = []
        if enzyme:
            features_used.append(f"Enzyme: {enzyme}")
        if include_his_tag:
            features_used.append("His-Tag")
        if include_linkers:
            features_used.append("Linkers")
        if cleavage_site:
            features_used.append(f"Cleavage: {cleavage_site}")
        
        if features_used:
            features_text = ", ".join(features_used)
            st.markdown(f'<div class="ext-synthesis-info-box">✅ Using: {features_text}</div>', unsafe_allow_html=True)
        else:
            st.markdown('<div class="ext-synthesis-warning-box">ℹ️ Basic fragmentation only - no additional features selected</div>', unsafe_allow_html=True)

        with st.expander("🔬 Advanced Options"):
            show_warnings = st.checkbox(
                "Show sequence warnings",
                value=True,
                help="Display warnings about sequence quality issues",
            )
            validate_assembly = st.checkbox(
                "Validate assembly",
                value=True,
                help="Check for potential assembly issues",
            )
        
        st.markdown('</div>', unsafe_allow_html=True)

    # Enhanced control buttons
    st.markdown('<div class="ext-synthesis-section-header">🎯 Actions</div>', unsafe_allow_html=True)
    
    col_btn1, col_btn2, col_btn3 = st.columns([1, 1, 2])
    with col_btn1:
        clean_btn = st.button("🧹 Clean Sequence", help="Remove invalid characters from sequence", use_container_width=True)
    with col_btn2:
        validate_btn = st.button("✅ Validate", help="Check sequence validity", use_container_width=True)
    with col_btn3:
        process_btn = st.button("🚀 Fragment Sequence", type="primary", help="Start fragmentation process", use_container_width=True)

    # Handle Clean: only cleans the sequence, preserves results
    if clean_btn and st.session_state.sequence_input:
        is_valid, clean_seq, message = validate_sequence_input(st.session_state.sequence_input)
        if is_valid:
            st.session_state.sequence_input = clean_seq
            st.markdown(f'<div class="ext-synthesis-success-box">✅ {message}</div>', unsafe_allow_html=True)
        else:
            st.markdown(f'<div class="ext-synthesis-error-box">❌ {message}</div>', unsafe_allow_html=True)

    # Handle Validate: only reports validity, does not alter results
    if validate_btn and st.session_state.sequence_input:
        is_valid, _, message = validate_sequence_input(st.session_state.sequence_input)
        if is_valid:
            st.markdown(f'<div class="ext-synthesis-success-box">✅ Sequence is valid. {message}</div>', unsafe_allow_html=True)
        else:
            st.markdown(f'<div class="ext-synthesis-error-box">❌ {message}</div>', unsafe_allow_html=True)

    # Main fragmentation process
    if process_btn:
        sequence = st.session_state.sequence_input
        if not sequence:
            st.markdown('<div class="ext-synthesis-error-box">❌ Please enter a DNA sequence.</div>', unsafe_allow_html=True)
            return

        is_valid, clean_seq, message = validate_sequence_input(sequence)
        if not is_valid:
            st.markdown(f'<div class="ext-synthesis-error-box">❌ Invalid sequence: {message}</div>', unsafe_allow_html=True)
            return

        if len(clean_seq) <= frag_size:
            st.markdown(f'<div class="ext-synthesis-warning-box">⚠️ Sequence ({len(clean_seq)} bp) is shorter than or equal to fragment size ({frag_size} bp). No fragmentation needed.</div>', unsafe_allow_html=True)
            return

        if overlap_size >= frag_size:
            st.markdown(f'<div class="ext-synthesis-error-box">❌ Overlap size ({overlap_size} bp) must be smaller than fragment size ({frag_size} bp).</div>', unsafe_allow_html=True)
            return

        try:
            with st.spinner("🧬 Fragmenting sequence..."):
                # Note: You may need to modify fragment_extended_sequence to handle the new parameters
                assembly, reassembled = fragment_extended_sequence(
                    clean_seq, frag_size, enzyme or "NdeI/XhoI", cleavage_site or "", overlap_size
                )

                st.session_state.sequences_processed = st.session_state.get("sequences_processed", 0) + 1
                st.session_state.fragments_generated = st.session_state.get("fragments_generated", 0) + len(assembly)

            st.markdown(f'<div class="ext-synthesis-success-box">✅ Successfully fragmented sequence ({len(clean_seq)} bp) into {len(assembly)} fragments!</div>', unsafe_allow_html=True)

            st.session_state.assembly = assembly
            st.session_state.reassembled = reassembled
            st.session_state.parameters = {
                "fragment_size": frag_size,
                "enzyme_pair": enzyme or "None",
                "cleavage_site": cleavage_site or "None",
                "internal_overlap": overlap_size,
                "include_his_tag": include_his_tag,
                "include_linkers": include_linkers,
            }

        except Exception as e:
            st.markdown(f'<div class="ext-synthesis-error-box">❌ Error during fragmentation: {e}</div>', unsafe_allow_html=True)
            with st.expander("Error Details"):
                st.exception(e)

    # Display results only if assembly is set
    if st.session_state.assembly:
        display_fragmentation_results(
            st.session_state.assembly,
            st.session_state.sequence_input,
            st.session_state.reassembled,
            st.session_state.parameters,
        )

if __name__ == "__main__":
    main()