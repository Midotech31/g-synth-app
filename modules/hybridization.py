# -*- coding: utf-8 -*-
"""
modules/hybridization.py
OLIGOSYNTHESIZER VERSION - F/R Hybridization Analysis

Features:
- Proper F (5'→3') and R (5'→3') hybridization analysis
- Watson-Crick complementarity checking
- Sticky end detection for direct ligation
- Interactive visualization showing true hybridization
- Modern blue theme interface
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
from collections import Counter

# ───────────────────────────────────────────────────────────────────────────────
# CORE FUNCTIONS FOR OLIGOSYNTHESIZER SEQUENCES
# ───────────────────────────────────────────────────────────────────────────────

def clean_dna_sequence(seq):
    """Clean and validate DNA sequence"""
    if not seq:
        return ""
    clean_seq = re.sub(r'\s+', '', seq.upper())
    clean_seq = re.sub(r'[^ATGC]', '', clean_seq)  # Only standard bases
    return clean_seq

def watson_crick_complement(base):
    """Get Watson-Crick complement of a base"""
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return complement_map.get(base, 'N')

def is_watson_crick_complement(base1, base2):
    """Check if two bases are Watson-Crick complementary"""
    return watson_crick_complement(base1) == base2

def analyze_fr_hybridization(f_seq, r_seq):
    """
    Analyze hybridization between F (5'→3') and R (5'→3') sequences
    F hybridizes with reverse of R following Watson-Crick rules
    """
    # R sequence reversed (3'→5' relative to F)
    r_reversed = r_seq[::-1]
    
    # Find optimal alignment for Watson-Crick complementarity
    best_alignment = find_optimal_alignment(f_seq, r_reversed)
    
    if not best_alignment:
        return None
    
    # Analyze the hybridization regions
    analysis = analyze_hybridization_regions(f_seq, r_reversed, best_alignment)
    
    # Add original sequences for reference
    analysis['original_sequences'] = {
        'f_seq': f_seq,
        'r_seq': r_seq,
        'r_reversed': r_reversed
    }
    
    return analysis

def find_optimal_alignment(f_seq, r_reversed):
    """Find optimal alignment for Watson-Crick hybridization"""
    best_score = 0
    best_alignment = None
    
    # Try different shifts to find best alignment
    max_shift = max(len(f_seq), len(r_reversed))
    
    for shift in range(-max_shift, max_shift + 1):
        score = 0
        matches = 0
        total_compared = 0
        
        # Compare sequences at this shift
        for i in range(len(f_seq)):
            j = i - shift
            if 0 <= j < len(r_reversed):
                total_compared += 1
                if is_watson_crick_complement(f_seq[i], r_reversed[j]):
                    score += 2
                    matches += 1
                else:
                    score -= 1  # Penalty for mismatch
        
        # Only consider alignments with reasonable overlap
        if total_compared >= min(len(f_seq), len(r_reversed)) * 0.3:
            if score > best_score:
                best_score = score
                best_alignment = {
                    'shift': shift,
                    'score': score,
                    'matches': matches,
                    'total_compared': total_compared,
                    'complementarity': (matches / total_compared * 100) if total_compared > 0 else 0
                }
    
    return best_alignment

def analyze_hybridization_regions(f_seq, r_reversed, alignment):
    """Analyze core and sticky end regions"""
    shift = alignment['shift']
    
    # Determine overlap region
    f_start = max(0, shift)
    f_end = min(len(f_seq), shift + len(r_reversed))
    
    r_start = max(0, -shift)
    r_end = min(len(r_reversed), len(f_seq) - shift)
    
    # Core hybridization region (where sequences overlap)
    f_core = f_seq[f_start:f_end]
    r_core = r_reversed[r_start:r_start + len(f_core)]
    
    # Calculate core complementarity
    core_matches = 0
    core_mismatches = 0
    core_details = []
    
    for i in range(len(f_core)):
        if i < len(r_core):
            f_base = f_core[i]
            r_base = r_core[i]
            is_complement = is_watson_crick_complement(f_base, r_base)
            
            core_details.append({
                'position': i,
                'f_base': f_base,
                'r_base': r_base,
                'complement': is_complement,
                'expected': watson_crick_complement(f_base)
            })
            
            if is_complement:
                core_matches += 1
            else:
                core_mismatches += 1
    
    core_complementarity = (core_matches / len(f_core) * 100) if len(f_core) > 0 else 0
    
    # Sticky end regions (overhangs)
    sticky_ends = {
        'f_5_prime': f_seq[:f_start] if f_start > 0 else "",
        'f_3_prime': f_seq[f_end:] if f_end < len(f_seq) else "",
        'r_5_prime': r_reversed[:r_start] if r_start > 0 else "",
        'r_3_prime': r_reversed[r_start + len(f_core):] if r_start + len(f_core) < len(r_reversed) else ""
    }
    
    analysis = {
        'core_region': {
            'f_core': f_core,
            'r_core': r_core,
            'length': len(f_core),
            'matches': core_matches,
            'mismatches': core_mismatches,
            'complementarity': core_complementarity,
            'details': core_details
        },
        'sticky_ends': sticky_ends,
        'alignment': alignment,
        'full_alignment': create_full_alignment(f_seq, r_reversed, shift)
    }
    
    return analysis

def create_full_alignment(f_seq, r_reversed, shift):
    """Create full alignment with gaps for visualization"""
    # Determine the full range
    start_pos = min(0, shift)
    end_pos = max(len(f_seq), shift + len(r_reversed))
    
    f_aligned = []
    r_aligned = []
    
    for pos in range(start_pos, end_pos):
        # Position in F sequence
        f_pos = pos
        # Position in R sequence
        r_pos = pos - shift
        
        if 0 <= f_pos < len(f_seq):
            f_aligned.append(f_seq[f_pos])
        else:
            f_aligned.append('-')
        
        if 0 <= r_pos < len(r_reversed):
            r_aligned.append(r_reversed[r_pos])
        else:
            r_aligned.append('-')
    
    return {
        'f_aligned': ''.join(f_aligned),
        'r_aligned': ''.join(r_aligned),
        'start_pos': start_pos
    }

def create_hybridization_visualization(analysis):
    """Create interactive hybridization visualization"""
    
    full_align = analysis['full_alignment']
    core = analysis['core_region']
    sticky = analysis['sticky_ends']
    
    f_aligned = full_align['f_aligned']
    r_aligned = full_align['r_aligned']
    
    positions = list(range(len(f_aligned)))
    
    # Color coding and hover text
    colors = []
    hover_texts = []
    
    # Determine region boundaries
    f_5_len = len(sticky['f_5_prime'])
    r_5_len = len(sticky['r_5_prime'])
    core_len = len(core['f_core'])
    
    for i, (f_base, r_base) in enumerate(zip(f_aligned, r_aligned)):
        if f_base == '-' or r_base == '-':
            # Sticky end region
            if f_base == '-':
                colors.append('#FF9800')  # Orange for R sticky end
                hover_texts.append(f"Position {i+1}: R sticky end ({r_base})")
            else:
                colors.append('#9C27B0')  # Purple for F sticky end
                hover_texts.append(f"Position {i+1}: F sticky end ({f_base})")
        else:
            # Core region - check complementarity
            if is_watson_crick_complement(f_base, r_base):
                colors.append('#2563EB')  # Blue for Watson-Crick complement
                hover_texts.append(f"Position {i+1}: {f_base}-{r_base} (Complementary)")
            else:
                colors.append('#FF6B6B')  # Red for mismatch
                expected = watson_crick_complement(f_base)
                hover_texts.append(f"Position {i+1}: {f_base}-{r_base} (Mismatch! Expected {f_base}-{expected})")
    
    # Create figure with 3 rows
    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=[
            'Forward Strand (5′→3′)', 
            'Reverse Strand (3′→5′)', 
            'Watson-Crick Complementarity'
        ],
        vertical_spacing=0.15
    )
    
    # Forward strand
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode='markers+text',
            marker=dict(size=14, color=colors, line=dict(width=1, color='white')),
            text=list(f_aligned),
            textposition='middle center',
            textfont=dict(size=10, color='white', family='monospace'),
            hovertext=hover_texts,
            hovertemplate='%{hovertext}<extra></extra>',
            name='Forward (5′→3′)',
            showlegend=False
        ),
        row=1, col=1
    )
    
    # Reverse strand (showing 3'→5' direction)
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode='markers+text',
            marker=dict(size=14, color=colors, line=dict(width=1, color='white')),
            text=list(r_aligned),
            textposition='middle center',
            textfont=dict(size=10, color='white', family='monospace'),
            hovertext=hover_texts,
            hovertemplate='%{hovertext}<extra></extra>',
            name='Reverse (3′→5′)',
            showlegend=False
        ),
        row=2, col=1
    )
    
    # Complementarity line graph
    comp_values = []
    for i, (f_base, r_base) in enumerate(zip(f_aligned, r_aligned)):
        if f_base == '-' or r_base == '-':
            comp_values.append(None)  # No value for sticky ends
        elif is_watson_crick_complement(f_base, r_base):
            comp_values.append(100)  # Perfect complement
        else:
            comp_values.append(0)  # Mismatch
    
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=comp_values,
            mode='lines+markers',
            line=dict(color='#2563EB', width=3),
            marker=dict(size=8, color='#1E40AF'),
            name='Complementarity',
            hovertemplate='Position %{x}: %{y}% Watson-Crick match<extra></extra>',
            connectgaps=False
        ),
        row=3, col=1
    )
    
    fig.update_layout(
        height=700,
        title="Oligosynthesizer F/R Hybridization Analysis",
        plot_bgcolor='#f8f9ff',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12)
    )
    
    # Update axes
    fig.update_yaxes(showticklabels=False, range=[0.5, 1.5], row=1, col=1)
    fig.update_yaxes(showticklabels=False, range=[0.5, 1.5], row=2, col=1)
    fig.update_yaxes(title_text="Complementarity (%)", range=[-10, 110], row=3, col=1)
    fig.update_xaxes(title_text="Position", row=3, col=1)
    
    return fig

# ───────────────────────────────────────────────────────────────────────────────
# MAIN APPLICATION
# ───────────────────────────────────────────────────────────────────────────────

def main():
    """Oligosynthesizer F/R hybridization analysis"""
    
    # Header
    st.markdown("""
    <div style="
        background: linear-gradient(135deg, #2563EB 0%, #1E40AF 100%);
        padding: 2rem;
        border-radius: 15px;
        text-align: center;
        margin-bottom: 2rem;
        color: white;
        box-shadow: 0 8px 32px rgba(37, 99, 235, 0.37);
    ">
        <h1 style="margin: 0; font-size: 2.5rem; font-weight: 800;">🧬 Oligosynthesizer F/R Hybridization</h1>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.2rem; opacity: 0.9;">
            Analyze Watson-Crick hybridization of F (5'→3') and R (5'→3') with sticky ends
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Initialize session state
    if "f_input" not in st.session_state:
        st.session_state["f_input"] = ""
    if "r_input" not in st.session_state:
        st.session_state["r_input"] = ""
    
    # Sidebar with explanation
    with st.sidebar:
        st.markdown("### 🧬 How It Works")
        st.info("""
        **Oligosynthesizer Sequences:**
        - F and R are both 5'→3'
        - F hybridizes with reverse of R
        - Watson-Crick complementarity (A-T, G-C)
        - Sticky ends included for ligation
        - Core region should be fully complementary
        """)
        
        st.markdown("### 📋 Examples")
        
        if st.button("Perfect Oligosynthesizer Pair"):
            st.session_state["f_input"] = "GATCATGCGTAAGCTTGATATCGAATTCCTGCA"
            st.session_state["r_input"] = "TGCAGGAATTCGATATCAAGCTTACGCATGATC"
            st.experimental_rerun()
        
        if st.button("With Sticky Ends"):
            st.session_state["f_input"] = "AATTGCGTAAGCTTGATATCGAATTCCTGC"
            st.session_state["r_input"] = "GCAGGAATTCGATATCAAGCTTACGCTTAA"
            st.experimental_rerun()
        
        if st.button("Clear All"):
            st.session_state["f_input"] = ""
            st.session_state["r_input"] = ""
            st.experimental_rerun()
    
    # Input section
    st.markdown("### 📝 Oligosynthesizer Sequence Input")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Forward Strand (F) - 5'→3'**")
        f_seq_raw = st.text_area(
            "",
            value=st.session_state["f_input"],
            key="f_input",
            height=150,
            placeholder="GATCATGCGTAAGCTTGATATCGAATTCCTGCA...",
            help="Enter Forward oligosynthesizer sequence (5'→3')"
        )
        
        if f_seq_raw:
            clean_f = clean_dna_sequence(f_seq_raw)
            gc_content = (clean_f.count('G') + clean_f.count('C')) / len(clean_f) * 100 if clean_f else 0
            st.info(f"📏 Length: {len(clean_f)} bp | 🧬 GC: {gc_content:.1f}%")
    
    with col2:
        st.markdown("**Reverse Strand (R) - 5'→3'**")
        r_seq_raw = st.text_area(
            "",
            value=st.session_state["r_input"],
            key="r_input",
            height=150,
            placeholder="TGCAGGAATTCGATATCAAGCTTACGCATGATC...",
            help="Enter Reverse oligosynthesizer sequence (5'→3')"
        )
        
        if r_seq_raw:
            clean_r = clean_dna_sequence(r_seq_raw)
            gc_content = (clean_r.count('G') + clean_r.count('C')) / len(clean_r) * 100 if clean_r else 0
            st.info(f"📏 Length: {len(clean_r)} bp | 🧬 GC: {gc_content:.1f}%")
    
    # Analysis button
    if st.button("🚀 Analyze F/R Hybridization", type="primary"):
        f_seq = clean_dna_sequence(f_seq_raw or "")
        r_seq = clean_dna_sequence(r_seq_raw or "")
        
        if not f_seq or not r_seq:
            st.error("⚠️ Please provide valid DNA sequences for both F and R strands.")
        else:
            with st.spinner("🔄 Analyzing oligosynthesizer F/R hybridization..."):
                
                # Analyze hybridization
                analysis = analyze_fr_hybridization(f_seq, r_seq)
                
                if not analysis:
                    st.error("❌ Could not find suitable hybridization alignment.")
                else:
                    st.success("✅ F/R hybridization analysis completed!")
                    
                    # Key metrics
                    core = analysis['core_region']
                    sticky = analysis['sticky_ends']
                    alignment = analysis['alignment']
                    
                    col_metric1, col_metric2, col_metric3, col_metric4 = st.columns(4)
                    
                    with col_metric1:
                        st.metric(
                            "Core Complementarity",
                            f"{core['complementarity']:.1f}%",
                            help="Watson-Crick complementarity in core region"
                        )
                    
                    with col_metric2:
                        st.metric(
                            "Core Length",
                            f"{core['length']} bp",
                            help="Length of hybridized core region"
                        )
                    
                    with col_metric3:
                        sticky_total = sum(len(s) for s in sticky.values())
                        st.metric(
                            "Total Sticky Ends",
                            f"{sticky_total} bp",
                            help="Total sticky end overhang length"
                        )
                    
                    with col_metric4:
                        st.metric(
                            "Overall Quality",
                            f"{alignment['complementarity']:.1f}%",
                            help="Overall alignment quality score"
                        )
                    
                    # Results tabs
                    tab1, tab2, tab3, tab4 = st.tabs([
                        "🎯 Interactive View", "🔍 Detailed Analysis", "📊 Results Table", "📋 Export"
                    ])
                    
                    with tab1:
                        st.markdown("#### Interactive F/R Hybridization Visualization")
                        
                        # Create and display visualization
                        viz_fig = create_hybridization_visualization(analysis)
                        st.plotly_chart(viz_fig, use_container_width=True)
                        
                        # Legend and explanation
                        st.markdown("#### Visualization Legend")
                        
                        col_leg1, col_leg2, col_leg3, col_leg4 = st.columns(4)
                        
                        with col_leg1:
                            st.markdown("🔵 **Watson-Crick Pairs**")
                            st.caption("Perfect A-T, G-C complementarity")
                        
                        with col_leg2:
                            st.markdown("🔴 **Mismatches**")
                            st.caption("Non-complementary base pairs")
                        
                        with col_leg3:
                            st.markdown("🟠 **R Sticky Ends**")
                            st.caption("R sequence overhangs")
                        
                        with col_leg4:
                            st.markdown("🟣 **F Sticky Ends**")
                            st.caption("F sequence overhangs")
                        
                        # Sequence orientation explanation
                        st.markdown("#### Sequence Orientation")
                        st.info("""
                        **F**: 5' → 3' (as synthesized)  
                        **R**: 3' → 5' (reversed for hybridization)  
                        **Core**: Watson-Crick complementary region  
                        **Sticky Ends**: Overhangs for direct ligation
                        """)
                    
                    with tab2:
                        st.markdown("#### Detailed Hybridization Analysis")
                        
                        # Core region assessment
                        if core['complementarity'] >= 98:
                            st.success(f"🟢 **Excellent Watson-Crick Complementarity** ({core['complementarity']:.1f}%)")
                            st.write("✅ Perfect oligosynthesizer design for hybridization")
                        elif core['complementarity'] >= 90:
                            st.warning(f"🟡 **Good Complementarity** ({core['complementarity']:.1f}%)")
                            st.write(f"⚠️ {core['mismatches']} mismatches in core region")
                        else:
                            st.error(f"🔴 **Poor Complementarity** ({core['complementarity']:.1f}%)")
                            st.write(f"❌ {core['mismatches']} mismatches - check oligosynthesizer design")
                        
                        # Core sequences
                        st.markdown("##### Core Hybridization Region")
                        
                        col_core1, col_core2 = st.columns(2)
                        
                        with col_core1:
                            st.markdown("**F Core (5'→3'):**")
                            st.code(core['f_core'], language=None)
                        
                        with col_core2:
                            st.markdown("**R Core (3'→5'):**")
                            st.code(core['r_core'], language=None)
                        
                        # Mismatch details
                        if core['mismatches'] > 0:
                            st.markdown("##### Mismatch Analysis")
                            
                            mismatch_data = []
                            for detail in core['details']:
                                if not detail['complement']:
                                    mismatch_data.append({
                                        'Position': detail['position'] + 1,
                                        'F Base': detail['f_base'],
                                        'R Base': detail['r_base'],
                                        'Expected R': detail['expected'],
                                        'Issue': f"Found {detail['r_base']}, expected {detail['expected']}"
                                    })
                            
                            if mismatch_data:
                                mismatch_df = pd.DataFrame(mismatch_data)
                                st.dataframe(mismatch_df, use_container_width=True)
                        
                        # Sticky ends analysis
                        st.markdown("##### Sticky End Analysis")
                        
                        sticky_data = []
                        if sticky['f_5_prime']:
                            sticky_data.append({
                                "Location": "F 5' overhang", 
                                "Sequence": sticky['f_5_prime'], 
                                "Length": len(sticky['f_5_prime'])
                            })
                        if sticky['f_3_prime']:
                            sticky_data.append({
                                "Location": "F 3' overhang", 
                                "Sequence": sticky['f_3_prime'], 
                                "Length": len(sticky['f_3_prime'])
                            })
                        if sticky['r_5_prime']:
                            sticky_data.append({
                                "Location": "R 5' overhang", 
                                "Sequence": sticky['r_5_prime'], 
                                "Length": len(sticky['r_5_prime'])
                            })
                        if sticky['r_3_prime']:
                            sticky_data.append({
                                "Location": "R 3' overhang", 
                                "Sequence": sticky['r_3_prime'], 
                                "Length": len(sticky['r_3_prime'])
                            })
                        
                        if sticky_data:
                            sticky_df = pd.DataFrame(sticky_data)
                            st.dataframe(sticky_df, use_container_width=True)
                            
                            if len(sticky_data) > 0:
                                st.success("✅ Sticky ends detected - ready for direct ligation after hybridization")
                        else:
                            st.info("ℹ️ No sticky ends detected - blunt end hybridization")
                    
                    with tab3:
                        st.markdown("#### Complete Analysis Results")
                        
                        # Comprehensive metrics table
                        orig_seqs = analysis['original_sequences']
                        
                        results_data = pd.DataFrame({
                            'Property': [
                                'F Sequence Length',
                                'R Sequence Length', 
                                'Core Hybridization Length',
                                'Core Watson-Crick Matches',
                                'Core Mismatches',
                                'Core Complementarity',
                                'Total Sticky End Length',
                                'Alignment Shift',
                                'Overall Quality Score'
                            ],
                            'Value': [
                                f"{len(orig_seqs['f_seq'])} bp",
                                f"{len(orig_seqs['r_seq'])} bp",
                                f"{core['length']} bp",
                                f"{core['matches']}",
                                f"{core['mismatches']}",
                                f"{core['complementarity']:.1f}%",
                                f"{sum(len(s) for s in sticky.values())} bp",
                                f"{alignment['shift']} bp",
                                f"{alignment['complementarity']:.1f}%"
                            ]
                        })
                        
                        st.dataframe(results_data, use_container_width=True)
                        
                        # Design quality assessment
                        st.markdown("#### Oligosynthesizer Design Quality")
                        
                        quality_score = core['complementarity']
                        
                        if quality_score >= 98:
                            st.success("🎉 **Excellent oligosynthesizer design** - Ready for synthesis and hybridization")
                        elif quality_score >= 90:
                            st.warning("⚠️ **Good design with minor issues** - Consider optimizing mismatched positions")
                        else:
                            st.error("🚨 **Design needs improvement** - Multiple mismatches detected")
                        
                        # Recommendations
                        recommendations = []
                        
                        if core['mismatches'] > 0:
                            recommendations.append(f"• Fix {core['mismatches']} mismatch(es) in core region for optimal hybridization")
                        
                        if core['length'] < min(len(orig_seqs['f_seq']), len(orig_seqs['r_seq'])) * 0.7:
                            recommendations.append("• Consider increasing core hybridization region")
                        
                        sticky_total = sum(len(s) for s in sticky.values())
                        if sticky_total == 0:
                            recommendations.append("• Add sticky ends if direct ligation is needed")
                        elif sticky_total > 20:
                            recommendations.append("• Consider shorter sticky ends for more efficient ligation")
                        
                        if recommendations:
                            st.markdown("##### Optimization Recommendations")
                            for rec in recommendations:
                                st.info(rec)
                        else:
                            st.success("✅ No optimization recommendations - excellent design!")
                    
                    with tab4:
                        st.markdown("#### Export Analysis Results")
                        
                        # Generate report
                        orig_seqs = analysis['original_sequences']
                        full_align = analysis['full_alignment']
                        
                        report = f"""
OLIGOSYNTHESIZER F/R HYBRIDIZATION ANALYSIS
==========================================

Input Sequences (5'→3'):
F: {orig_seqs['f_seq']}
R: {orig_seqs['r_seq']}

Hybridization Analysis:
F (5'→3'): {full_align['f_aligned']}
R (3'→5'): {full_align['r_aligned']}

Core Region Analysis:
- Core F: {core['f_core']}
- Core R: {core['r_core']}
- Length: {core['length']} bp
- Watson-Crick Matches: {core['matches']}
- Mismatches: {core['mismatches']}
- Complementarity: {core['complementarity']:.1f}%

Sticky End Analysis:
- F 5' overhang: {sticky['f_5_prime'] or 'None'}
- F 3' overhang: {sticky['f_3_prime'] or 'None'}
- R 5' overhang: {sticky['r_5_prime'] or 'None'}
- R 3' overhang: {sticky['r_3_prime'] or 'None'}
- Total sticky length: {sum(len(s) for s in sticky.values())} bp

Quality Assessment:
Core Complementarity: {core['complementarity']:.1f}%
Overall Score: {alignment['complementarity']:.1f}%
Design Quality: {"Excellent" if core['complementarity'] >= 98 else "Good" if core['complementarity'] >= 90 else "Needs Improvement"}

Generated by Oligosynthesizer F/R Hybridization Analyzer
"""
                        
                        st.text_area("Complete Analysis Report", report, height=400)
                        
                        # Download options
                        col_dl1, col_dl2, col_dl3 = st.columns(3)
                        
                        with col_dl1:
                            st.download_button(
                                "📄 Download Report",
                                report,
                                "fr_hybridization_report.txt",
                                "text/plain"
                            )
                        
                        with col_dl2:
                            core_fasta = f">F_Core_5to3\n{core['f_core']}\n>R_Core_3to5\n{core['r_core']}"
                            st.download_button(
                                "📄 Download Core FASTA",
                                core_fasta,
                                "core_hybridization.fasta",
                                "text/plain"
                            )
                        
                        with col_dl3:
                            full_fasta = f">F_Original_5to3\n{orig_seqs['f_seq']}\n>R_Original_5to3\n{orig_seqs['r_seq']}"
                            st.download_button(
                                "📄 Download Original FASTA",
                                full_fasta,
                                "original_fr_sequences.fasta",
                                "text/plain"
                            )

def app():
    """Alternative entrypoint function for compatibility"""
    main()

if __name__ == "__main__":
    main()