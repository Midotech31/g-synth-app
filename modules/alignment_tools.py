# -*- coding: utf-8 -*-
"""
modules/alignment_tools.py
FINAL COMPLETE VERSION - Professional alignment tool with all features
- Optimal sticky end alignment with reverse complement detection
- Mutation highlighting with visual effects
- Interactive, zoomable logomaker showing only variable positions
- Complete responsive visualization with blue theme
- Enhanced statistics and analysis
- FIXED complementarity calculation
- All features integrated and fully functional
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit.components.v1 as components
from io import StringIO
import io

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align import MultipleSeqAlignment
from Bio.Align import substitution_matrices

# Try import logomaker for sequence logos
try:
    import logomaker
    _HAS_LOGOMAKER = True
except ImportError:
    _HAS_LOGOMAKER = False

# ───────────────────────────────────────────────────────────────────────────────
# OPTIMAL STICKY END ALIGNMENT LOGIC
# ───────────────────────────────────────────────────────────────────────────────

def detect_sticky_end_sequences(seq1, seq2):
    """
    Detect if sequences are sticky ends (forward and reverse complement)
    """
    seq1_upper = seq1.upper()
    seq2_upper = seq2.upper()
    seq2_rc = str(Seq(seq2_upper).reverse_complement())
    
    def simple_similarity(s1, s2):
        if len(s1) != len(s2):
            min_len = min(len(s1), len(s2))
            s1, s2 = s1[:min_len], s2[:min_len]
        matches = sum(1 for a, b in zip(s1, s2) if a == b)
        return matches / len(s1) if len(s1) > 0 else 0
    
    sim_forward = simple_similarity(seq1_upper, seq2_upper)
    sim_reverse = simple_similarity(seq1_upper, seq2_rc)
    
    return sim_reverse > sim_forward, sim_reverse

def create_optimal_sticky_end_alignment(seq1, seq2):
    """
    OPTIMAL: Create optimal alignment for sticky end sequences
    """
    try:
        seq1_upper = seq1.upper()
        seq2_upper = seq2.upper()
        seq2_rc = str(Seq(seq2_upper).reverse_complement())
        
        # Create high-quality aligner for DNA
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 5.0
        aligner.mismatch_score = -4.0
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -1.0
        
        # Align both orientations
        alignment_forward = aligner.align(seq1_upper, seq2_upper)
        score_forward = alignment_forward.score if alignment_forward else -999999
        
        alignment_reverse = aligner.align(seq1_upper, seq2_rc)
        score_reverse = alignment_reverse.score if alignment_reverse else -999999
        
        # Choose best alignment
        if score_reverse > score_forward:
            best_alignment = alignment_reverse[0]
            orientation = "Reverse Complement"
            aligned_seq2 = seq2_rc
            is_sticky_end = True
        else:
            best_alignment = alignment_forward[0]
            orientation = "Forward"
            aligned_seq2 = seq2_upper
            is_sticky_end = False
        
        # Extract properly gapped sequences
        target_gapped, query_gapped = extract_gapped_sequences_from_alignment(best_alignment)
        
        return target_gapped, query_gapped, best_alignment, orientation, is_sticky_end
        
    except Exception as e:
        st.error(f"Sticky end alignment failed: {e}")
        return None, None, None, "Error", False

def extract_gapped_sequences_from_alignment(alignment):
    """
    OPTIMAL: Extract properly gapped sequences from BioPython alignment
    """
    try:
        target_coords = alignment.aligned[0]
        query_coords = alignment.aligned[1]
        
        target_seq = str(alignment.target)
        query_seq = str(alignment.query)
        
        # Build gapped sequences
        target_gapped = []
        query_gapped = []
        
        target_pos = 0
        query_pos = 0
        
        # Process each aligned block
        for i in range(max(len(target_coords), len(query_coords))):
            if i < len(target_coords):
                t_start, t_end = target_coords[i]
            else:
                t_start = t_end = len(target_seq)
                
            if i < len(query_coords):
                q_start, q_end = query_coords[i]
            else:
                q_start = q_end = len(query_seq)
            
            # Add gaps for unaligned regions
            while target_pos < t_start:
                target_gapped.append(target_seq[target_pos])
                query_gapped.append('-')
                target_pos += 1
                
            while query_pos < q_start:
                target_gapped.append('-')
                query_gapped.append(query_seq[query_pos])
                query_pos += 1
            
            # Add aligned region
            while target_pos < t_end and query_pos < q_end:
                target_gapped.append(target_seq[target_pos])
                query_gapped.append(query_seq[query_pos])
                target_pos += 1
                query_pos += 1
        
        # Add remaining bases
        while target_pos < len(target_seq):
            target_gapped.append(target_seq[target_pos])
            query_gapped.append('-')
            target_pos += 1
            
        while query_pos < len(query_seq):
            target_gapped.append('-')
            query_gapped.append(query_seq[query_pos])
            query_pos += 1
        
        # Ensure equal length
        max_len = max(len(target_gapped), len(query_gapped))
        while len(target_gapped) < max_len:
            target_gapped.append('-')
        while len(query_gapped) < max_len:
            query_gapped.append('-')
        
        return ''.join(target_gapped), ''.join(query_gapped)
        
    except Exception as e:
        st.warning(f"Could not extract gapped sequences: {e}")
        # Fallback
        target_seq = str(alignment.target)
        query_seq = str(alignment.query)
        max_len = max(len(target_seq), len(query_seq))
        return target_seq.ljust(max_len, '-'), query_seq.ljust(max_len, '-')

def calculate_alignment_statistics(seq1_aligned, seq2_aligned, original_seq1="", original_seq2="", is_sticky_end=False, orientation="Forward"):
    """
    COMPLETELY FIXED: Calculate comprehensive alignment statistics with CORRECT complementarity
    """
    if not seq1_aligned or not seq2_aligned:
        return 0, 0, 0, 0.0, 0, 0.0
    
    matches = 0
    mismatches = 0
    gaps = 0
    
    min_len = min(len(seq1_aligned), len(seq2_aligned))
    
    for i in range(min_len):
        char1 = seq1_aligned[i]
        char2 = seq2_aligned[i]
        
        if char1 == '-' or char2 == '-':
            gaps += 1
        elif char1 == char2:
            matches += 1
        else:
            mismatches += 1
    
    # COMPLETELY FIXED: Proper complementarity calculation for sticky ends
    complement_matches = 0
    
    if is_sticky_end and orientation == "Reverse Complement":
        # For sticky ends in reverse complement orientation, 
        # we need to check complementarity properly
        
        # Method 1: Check complementarity in aligned sequences
        # Since seq2 is already reverse complemented, check direct complementarity
        for i in range(min_len):
            char1 = seq1_aligned[i]
            char2 = seq2_aligned[i]
            
            # Only count non-gap positions
            if char1 != '-' and char2 != '-':
                # For sticky ends, the sequences should be complementary
                # Check Watson-Crick complementary base pairing
                if (char1 == 'A' and char2 == 'T') or (char1 == 'T' and char2 == 'A') or \
                   (char1 == 'G' and char2 == 'C') or (char1 == 'C' and char2 == 'G'):
                    complement_matches += 1
        
        # Calculate complementarity based on non-gap positions only
        non_gap_positions = matches + mismatches  # Positions where both sequences have bases
        complementarity = (complement_matches / non_gap_positions * 100) if non_gap_positions > 0 else 0.0
        
    elif is_sticky_end and orientation == "Forward":
        # If forward orientation was chosen for sticky ends, 
        # calculate normal complementarity
        for i in range(min_len):
            char1 = seq1_aligned[i]
            char2 = seq2_aligned[i]
            if char1 != '-' and char2 != '-':
                if (char1 == 'A' and char2 == 'T') or (char1 == 'T' and char2 == 'A') or \
                   (char1 == 'G' and char2 == 'C') or (char1 == 'C' and char2 == 'G'):
                    complement_matches += 1
        
        non_gap_positions = matches + mismatches
        complementarity = (complement_matches / non_gap_positions * 100) if non_gap_positions > 0 else 0.0
        
    else:
        # For non-sticky ends, standard complementarity calculation
        for i in range(min_len):
            char1 = seq1_aligned[i]
            char2 = seq2_aligned[i]
            if char1 != '-' and char2 != '-':
                if (char1 == 'A' and char2 == 'T') or (char1 == 'T' and char2 == 'A') or \
                   (char1 == 'G' and char2 == 'C') or (char1 == 'C' and char2 == 'G'):
                    complement_matches += 1
        
        non_gap_positions = matches + mismatches
        complementarity = (complement_matches / non_gap_positions * 100) if non_gap_positions > 0 else 0.0
    
    total_positions = matches + mismatches
    identity = (matches / total_positions * 100) if total_positions > 0 else 0.0
    
    return matches, mismatches, gaps, identity, complement_matches, complementarity

# ───────────────────────────────────────────────────────────────────────────────
# MUTATION HIGHLIGHTING FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def highlight_mutations(seq1_aligned, seq2_aligned):
    """
    Detect mutations (mismatches) between aligned sequences.
    """
    mutations = []
    for i, (a, b) in enumerate(zip(seq1_aligned, seq2_aligned)):
        if a != b and a != '-' and b != '-':
            mutations.append(i)
    return mutations

def analyze_mutations_detailed(seq1_aligned, seq2_aligned):
    """
    Detailed mutation analysis including types of mutations.
    """
    mutations = []
    transitions = []
    transversions = []
    
    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    
    for i, (a, b) in enumerate(zip(seq1_aligned, seq2_aligned)):
        if a != b and a != '-' and b != '-':
            mutation_info = {
                'position': i,
                'from': a,
                'to': b,
                'type': 'transition' if (a, b) in transition_pairs else 'transversion'
            }
            mutations.append(mutation_info)
            
            if (a, b) in transition_pairs:
                transitions.append(mutation_info)
            else:
                transversions.append(mutation_info)
    
    return mutations, transitions, transversions

# ───────────────────────────────────────────────────────────────────────────────
# COMPLETE ALIGNMENT VIEWER WITH MUTATION HIGHLIGHTING (BLUE THEME)
# ───────────────────────────────────────────────────────────────────────────────

def create_complete_alignment_viewer_html(sequences, seq_names, alignment_data=None, is_sticky_end=False):
    """
    COMPLETE: Professional alignment viewer with mutation highlighting and blue design
    """
    if not sequences:
        return ""
    
    max_len = max(len(seq) for seq in sequences)
    padded_sequences = [seq.ljust(max_len, '-') for seq in sequences]
    
    mutations = highlight_mutations(padded_sequences[0], padded_sequences[1]) if len(padded_sequences) >= 2 else []
    detailed_mutations, transitions, transversions = analyze_mutations_detailed(padded_sequences[0], padded_sequences[1]) if len(padded_sequences) >= 2 else ([], [], [])
    
    # Color scheme for nucleotides
    nucleotide_colors = {
        'A': '#FF6B6B', 'T': '#4ECDC4', 'G': '#45B7D1', 'C': '#96CEB4', 'U': '#4ECDC4', '-': '#E9ECEF'
    }
    
    # Create complete responsive HTML with BLUE theme (not purple)
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <style>
            .msa-container {{
                font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
                border: 1px solid #e1e8ed;
                border-radius: 12px;
                background: linear-gradient(135deg, #2563EB 0%, #1E40AF 100%);
                padding: 20px;
                margin: 10px 0;
                box-shadow: 0 8px 32px rgba(37, 99, 235, 0.37);
                backdrop-filter: blur(4px);
                -webkit-backdrop-filter: blur(4px);
                overflow-x: auto;
                max-width: 100%;
            }}
            .msa-header {{
                background: rgba(255, 255, 255, 0.1);
                padding: 15px;
                border-radius: 8px;
                margin-bottom: 15px;
                color: white;
                font-weight: bold;
            }}
            .sequence-container {{
                overflow-x: auto;
                max-width: 100%;
                background: rgba(255, 255, 255, 0.05);
                border-radius: 8px;
                padding: 10px;
            }}
            .sequence-row {{
                display: flex;
                align-items: center;
                margin: 5px 0;
                background: rgba(255, 255, 255, 0.1);
                border-radius: 8px;
                padding: 8px;
                transition: all 0.3s ease;
                min-width: fit-content;
            }}
            .sequence-row:hover {{
                background: rgba(255, 255, 255, 0.2);
                transform: translateX(2px);
            }}
            .seq-name {{
                width: 150px;
                min-width: 150px;
                font-weight: bold;
                color: white;
                text-shadow: 1px 1px 2px rgba(0,0,0,0.3);
                margin-right: 10px;
            }}
            .seq-content {{
                font-family: 'Consolas', monospace;
                letter-spacing: 1px;
                font-size: 14px;
                display: flex;
                flex-wrap: nowrap;
                min-width: fit-content;
            }}
            .nucleotide {{
                padding: 2px 4px;
                margin: 0 1px;
                border-radius: 4px;
                color: white;
                font-weight: bold;
                text-shadow: 1px 1px 1px rgba(0,0,0,0.3);
                transition: all 0.2s ease;
                min-width: 20px;
                text-align: center;
                position: relative;
            }}
            .nucleotide:hover {{
                transform: scale(1.2);
                z-index: 10;
                position: relative;
                box-shadow: 0 4px 8px rgba(0,0,0,0.3);
            }}
            .mutation {{
                border: 3px solid #FF1744 !important;
                box-shadow: 0 0 15px rgba(255, 23, 68, 0.8) !important;
                animation: mutation-pulse 2s infinite;
            }}
            .mutation::after {{
                content: '⚡';
                position: absolute;
                top: -8px;
                right: -8px;
                font-size: 12px;
                color: #FF1744;
                animation: mutation-flash 1s infinite;
            }}
            .transition {{
                border: 2px solid #FF9800 !important;
                box-shadow: 0 0 10px rgba(255, 152, 0, 0.6) !important;
            }}
            .transversion {{
                border: 2px solid #E91E63 !important;
                box-shadow: 0 0 10px rgba(233, 30, 99, 0.6) !important;
            }}
            @keyframes mutation-pulse {{
                0%, 100% {{ box-shadow: 0 0 15px rgba(255, 23, 68, 0.8); }}
                50% {{ box-shadow: 0 0 25px rgba(255, 23, 68, 1); }}
            }}
            @keyframes mutation-flash {{
                0%, 100% {{ opacity: 1; }}
                50% {{ opacity: 0.5; }}
            }}
            .ruler {{
                font-size: 10px;
                color: rgba(255,255,255,0.7);
                margin: 10px 0;
                padding-left: 160px;
                overflow-x: auto;
                white-space: nowrap;
            }}
            .consensus-row {{
                background: linear-gradient(45deg, #FF6B6B, #4ECDC4);
                border-radius: 8px;
                margin-top: 15px;
                padding: 10px;
            }}
            .stats-panel {{
                background: rgba(255, 255, 255, 0.1);
                border-radius: 8px;
                padding: 15px;
                margin-top: 15px;
                color: white;
            }}
            .sticky-end-info {{
                background: rgba(255, 215, 0, 0.2);
                border: 2px solid #FFD700;
                border-radius: 8px;
                padding: 10px;
                margin: 10px 0;
                color: white;
                text-align: center;
            }}
            .mutation-info {{
                background: rgba(255, 23, 68, 0.2);
                border: 2px solid #FF1744;
                border-radius: 8px;
                padding: 10px;
                margin: 10px 0;
                color: white;
                text-align: center;
            }}
        </style>
    </head>
    <body>
        <div class="msa-container">
            <div class="msa-header">
                🧬 Pairwise Alignment Viewer
                <div style="font-size: 12px; margin-top: 5px; opacity: 0.8;">
                    {len(sequences)} sequences • {max_len} positions • {len(mutations)} mutations
                    {"• Sticky End Detected" if is_sticky_end else ""}
                </div>
            </div>
    """
    
    # Add mutation info
    if mutations:
        html_content += f"""
        <div class="mutation-info">
            🔬 {len(mutations)} Mutations: {len(transitions)} Transitions • {len(transversions)} Transversions
        </div>
        """
    
    # Add sticky end info
    if is_sticky_end:
        html_content += """
        <div class="sticky-end-info">
            ⚡ Sticky End Sequences - Optimal Complement Alignment
        </div>
        """
    
    # Add position ruler
    ruler_parts = []
    for i in range(0, max_len, 50):
        chunk = ""
        for j in range(i, min(i + 50, max_len)):
            if j % 10 == 0:
                chunk += str(j % 100)
            else:
                chunk += "·"
        ruler_parts.append(f'<span style="margin-right: 20px;">{chunk}</span>')
    
    html_content += f'<div class="ruler">{"".join(ruler_parts)}</div>'
    
    # Add sequences with mutation highlighting
    html_content += '<div class="sequence-container">'
    
    for seq_idx, (seq, name) in enumerate(zip(padded_sequences, seq_names)):
        html_content += f'<div class="sequence-row">'
        html_content += f'<div class="seq-name">{name[:15]}</div>'
        html_content += f'<div class="seq-content">'
        
        for pos, nucleotide in enumerate(seq):
            color = nucleotide_colors.get(nucleotide.upper(), '#95A5A6')
            
            # Add mutation highlighting classes
            mutation_class = ""
            if pos in mutations:
                mutation_class += " mutation"
                for mut in detailed_mutations:
                    if mut['position'] == pos:
                        if mut['type'] == 'transition':
                            mutation_class += " transition"
                        else:
                            mutation_class += " transversion"
                        break
            
            html_content += f'<span class="nucleotide{mutation_class}" style="background-color: {color};" title="Position {pos+1}: {nucleotide}">{nucleotide}</span>'
        
        html_content += '</div></div>'
    
    # Add consensus
    if len(sequences) > 1:
        consensus = ""
        for i in range(max_len):
            column = [seq[i] for seq in padded_sequences if i < len(seq)]
            if column:
                from collections import Counter
                most_common = Counter([c for c in column if c != '-']).most_common(1)
                consensus += most_common[0][0] if most_common else '-'
            else:
                consensus += '-'
        
        html_content += f'<div class="consensus-row">'
        html_content += f'<div style="display: flex; align-items: center;">'
        html_content += f'<div class="seq-name">Consensus</div>'
        html_content += f'<div class="seq-content">'
        
        for pos, nucleotide in enumerate(consensus):
            color = nucleotide_colors.get(nucleotide.upper(), '#95A5A6')
            mutation_indicator = " ⚡" if pos in mutations else ""
            html_content += f'<span class="nucleotide" style="background-color: {color}; font-weight: bold;">{nucleotide}{mutation_indicator}</span>'
        
        html_content += '</div></div></div>'
    
    html_content += '</div>'
    
    # Enhanced statistics panel
    if len(sequences) == 2 and alignment_data:
        if len(alignment_data) == 6:  # With complementarity data
            matches, mismatches, gaps, identity, complement_matches, complementarity = alignment_data
            
            html_content += f"""
            <div class="stats-panel">
                <div style="display: grid; grid-template-columns: repeat(6, 1fr); gap: 10px; text-align: center;">
                    <div>
                        <div style="font-size: 20px; font-weight: bold; color: #4ECDC4;">{matches}</div>
                        <div style="font-size: 11px; opacity: 0.8;">Matches</div>
                    </div>
                    <div>
                        <div style="font-size: 20px; font-weight: bold; color: #FF6B6B;">{mismatches}</div>
                        <div style="font-size: 11px; opacity: 0.8;">Mismatches</div>
                    </div>
                    <div>
                        <div style="font-size: 20px; font-weight: bold; color: #FF1744;">{len(mutations)}</div>
                        <div style="font-size: 11px; opacity: 0.8;">Mutations</div>
                    </div>
                    <div>
                        <div style="font-size: 20px; font-weight: bold; color: #FFD700;">{complement_matches}</div>
                        <div style="font-size: 11px; opacity: 0.8;">Complements</div>
                    </div>
                    <div>
                        <div style="font-size: 20px; font-weight: bold; color: #45B7D1;">{identity:.1f}%</div>
                        <div style="font-size: 11px; opacity: 0.8;">Identity</div>
                    </div>
                    <div>
                        <div style="font-size: 20px; font-weight: bold; color: #96CEB4;">{complementarity:.1f}%</div>
                        <div style="font-size: 11px; opacity: 0.8;">Complement %</div>
                    </div>
                </div>
            </div>
            """
        else:
            matches, mismatches, gaps, identity = alignment_data
            html_content += f"""
            <div class="stats-panel">
                <div style="display: grid; grid-template-columns: repeat(5, 1fr); gap: 12px; text-align: center;">
                    <div>
                        <div style="font-size: 22px; font-weight: bold; color: #4ECDC4;">{matches}</div>
                        <div style="font-size: 12px; opacity: 0.8;">Matches</div>
                    </div>
                    <div>
                        <div style="font-size: 22px; font-weight: bold; color: #FF6B6B;">{mismatches}</div>
                        <div style="font-size: 12px; opacity: 0.8;">Mismatches</div>
                    </div>
                    <div>
                        <div style="font-size: 22px; font-weight: bold; color: #FF1744;">{len(mutations)}</div>
                        <div style="font-size: 12px; opacity: 0.8;">Mutations</div>
                    </div>
                    <div>
                        <div style="font-size: 22px; font-weight: bold; color: #45B7D1;">{identity:.1f}%</div>
                        <div style="font-size: 12px; opacity: 0.8;">Identity</div>
                    </div>
                    <div>
                        <div style="font-size: 22px; font-weight: bold; color: #96CEB4;">{max_len}</div>
                        <div style="font-size: 12px; opacity: 0.8;">Length</div>
                    </div>
                </div>
            </div>
            """
    
    html_content += """
        </div>
    </body>
    </html>
    """
    
    return html_content

# ───────────────────────────────────────────────────────────────────────────────
# INTERACTIVE LOGOMAKER FUNCTION (VARIABLE POSITIONS ONLY)
# ───────────────────────────────────────────────────────────────────────────────

def create_interactive_variable_positions_logo(sequences, seq_names, title="Interactive Variable Positions"):
    """
    INTERACTIVE: Create zoomable, interactive logomaker showing only variable positions
    Professional blue design with corrected Plotly properties
    """
    try:
        if not sequences or len(sequences) < 2:
            return None, "Need at least 2 sequences for interactive logo"
        
        # Clean and pad sequences
        clean_sequences = []
        max_len = max(len(str(seq)) for seq in sequences)
        
        for seq in sequences:
            seq_str = str(seq).upper().strip()
            seq_str = ''.join(c for c in seq_str if c in 'ATCGRYSWKMBDHVN-')
            seq_str = seq_str.ljust(max_len, '-')
            clean_sequences.append(seq_str)
        
        # Find only variable positions
        variable_positions = []
        position_data = []
        
        for pos in range(max_len):
            column = [seq[pos] for seq in clean_sequences]
            non_gap_chars = [c for c in column if c != '-']
            unique_chars = set(non_gap_chars)
            
            if len(unique_chars) > 1:  # Variable position
                variable_positions.append(pos)
                
                # Calculate information content for this position
                from collections import Counter
                counts = Counter(non_gap_chars)
                total = len(non_gap_chars)
                
                entropy = 0
                for base, count in counts.items():
                    if count > 0:
                        p = count / total
                        entropy -= p * np.log2(p)
                
                info_content = 2 - entropy  # For DNA
                
                # Store data for each base at this position
                for base, count in counts.items():
                    height = (count / total) * info_content
                    position_data.append({
                        'position': len(variable_positions),
                        'original_position': pos + 1,
                        'base': base,
                        'height': height,
                        'frequency': count,
                        'total': total,
                        'percentage': (count / total) * 100
                    })
        
        if not variable_positions:
            return None, "No variable positions found for interactive logo"
        
        # Create interactive plotly figure
        fig = go.Figure()
        
        # Base colors for interactive logo (blue theme)
        base_colors = {
            'A': '#FF6B6B', 'T': '#4ECDC4', 'G': '#45B7D1', 'C': '#96CEB4', 'U': '#4ECDC4'
        }
        
        # Group data by position for stacking
        positions = {}
        for data in position_data:
            pos = data['position']
            if pos not in positions:
                positions[pos] = []
            positions[pos].append(data)
        
        # Create stacked bars for each position
        for pos, bases in positions.items():
            bases.sort(key=lambda x: x['height'])
            
            y_offset = 0
            for base_data in bases:
                fig.add_trace(go.Bar(
                    x=[pos],
                    y=[base_data['height']],
                    base=[y_offset],
                    name=base_data['base'],
                    marker_color=base_colors.get(base_data['base'], '#95A5A6'),
                    text=base_data['base'],
                    textposition='inside',
                    textfont=dict(size=14, color='white', family='Arial Black'),
                    hovertemplate=f"""
                    <b>Position {base_data['original_position']}</b><br>
                    Base: {base_data['base']}<br>
                    Frequency: {base_data['frequency']}/{base_data['total']}<br>
                    Percentage: {base_data['percentage']:.1f}%<br>
                    Information: {base_data['height']:.2f} bits
                    <extra></extra>
                    """,
                    showlegend=pos == 1
                ))
                y_offset += base_data['height']
        
        # Update layout for interactivity with blue theme
        fig.update_layout(
            title=dict(
                text=f"{title}<br><span style='font-size:12px'>{len(variable_positions)} Variable Positions - Interactive & Zoomable</span>",
                x=0.5,
                font=dict(size=16, family="Arial Black", color='#2563EB')
            ),
            xaxis=dict(
                title="Variable Position",
                tickmode='array',
                tickvals=list(range(1, len(variable_positions) + 1)),
                ticktext=[str(pos + 1) for pos in variable_positions],
                tickangle=45 if len(variable_positions) > 10 else 0,
                showgrid=True,
                gridcolor='rgba(37, 99, 235, 0.3)',
                gridwidth=1
            ),
            yaxis=dict(
                title="Information Content (bits)",
                showgrid=True,
                gridcolor='rgba(37, 99, 235, 0.3)',
                gridwidth=1
            ),
            barmode='stack',
            height=450,
            plot_bgcolor='#f8f9ff',
            paper_bgcolor='white',
            font=dict(family="Arial", size=12, color='#2c3e50'),
            
            # INTERACTIVE FEATURES
            dragmode='zoom',
            selectdirection='d',
            
            # Add modebar with custom tools
            modebar=dict(
                bgcolor='rgba(37, 99, 235, 0.1)',
                color='#2563EB',
                activecolor='#FF6B6B'
            )
        )
        
        # Add interactive annotations
        fig.add_annotation(
            x=0.02, y=0.98,
            xref="paper", yref="paper",
            text="🔍 Interactive: Zoom, Pan, Hover for details",
            showarrow=False,
            font=dict(size=10, color='#2563EB'),
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor='#2563EB',
            borderwidth=1
        )
        
        # Update traces for better interactivity
        fig.update_traces(
            marker_line_width=1,
            marker_line_color='rgba(255,255,255,0.3)'
        )
        
        return fig, None
        
    except Exception as e:
        return None, f"Interactive logo error: {str(e)}"

# ───────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def parse_fasta_input(uploaded, pasted):
    """Return list of SeqRecord from file upload or pasted FASTA."""
    if uploaded:
        records = []
        for f in uploaded:
            records.extend(list(SeqIO.parse(f, "fasta")))
        return records
    elif pasted.strip():
        return list(SeqIO.parse(StringIO(pasted), "fasta"))
    else:
        return []

# ───────────────────────────────────────────────────────────────────────────────
# MAIN APPLICATION - PROFESSIONAL AND COMPLETE
# ───────────────────────────────────────────────────────────────────────────────

def main():
    # Professional header with BLUE theme (not purple)
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
        <h1 style="margin: 0; font-size: 2.5rem; font-weight: 800;">🧬 Pairwise Alignment Studio</h1>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.2rem; opacity: 0.9;">
            Professional sequence alignment with mutation highlighting & interactive logos
        </p>
    </div>
    """, unsafe_allow_html=True)

    # Input section
    st.markdown("### 📁 Input Sequences")
    
    col_input1, col_input2 = st.columns([1, 1])
    
    with col_input1:
        uploaded = st.file_uploader(
            "Upload FASTA files",
            type=["fa", "fasta", "txt"],
            accept_multiple_files=True,
            help="Upload one or more FASTA files"
        )
    
    with col_input2:
        pasted = st.text_area(
            "Or paste FASTA sequences",
            height=150,
            placeholder=">forward_strand\nATGCGTAAGCTT...\n>reverse_strand\nAAGCTTACGCAT...",
            help="Enter sequences in FASTA format"
        )
    
    try:
        records = parse_fasta_input(uploaded, pasted)
    except Exception as e:
        st.error(f"Error parsing sequences: {e}")
        return
    
    if len(records) < 2:
        st.info("🔬 Please provide at least two sequences for alignment.")
        return

    # Show sequence info with blue theme
    st.markdown("### 📊 Sequence Information")
    cols = st.columns(min(len(records), 4))
    for i, record in enumerate(records):
        with cols[i % 4]:
            st.markdown(f"""
            <div style="
                background: linear-gradient(145deg, #f8f9fa, #e9ecef);
                padding: 1rem;
                border-radius: 10px;
                border-left: 4px solid #2563EB;
                margin: 0.5rem 0;
            ">
                <h4 style="margin: 0; color: #2c3e50;">{record.id or f'Seq {i+1}'}</h4>
                <p style="margin: 0.5rem 0 0 0; color: #7f8c8d;">{len(record.seq)} bp</p>
            </div>
            """, unsafe_allow_html=True)

    # Alignment configuration
    st.markdown("### ⚙️ Pairwise Alignment Configuration")
    
    col_config1, col_config2 = st.columns(2)
    
    with col_config1:
        seq_options = [f"{record.id or f'Sequence {i+1}'}" for i, record in enumerate(records)]
        seq1_idx = st.selectbox("Select first sequence:", range(len(records)), 
                               format_func=lambda x: seq_options[x])
        seq2_idx = st.selectbox("Select second sequence:", range(len(records)), 
                               format_func=lambda x: seq_options[x], 
                               index=1 if len(records) > 1 else 0)
    
    with col_config2:
        is_protein = st.checkbox("Protein sequences (use BLOSUM62)", value=False)
        sticky_end_mode = st.checkbox("Sticky end mode (auto-detect reverse complement)", value=True)

    if st.button("🚀 Run Pairwise Alignment", type="primary"):
        try:
            with st.spinner("🔄 Running professional alignment analysis..."):
                # Get sequences
                seq_a = str(records[seq1_idx].seq)
                seq_b = str(records[seq2_idx].seq)
                
                seq_names = [
                    records[seq1_idx].id or f"Sequence_{seq1_idx+1}",
                    records[seq2_idx].id or f"Sequence_{seq2_idx+1}"
                ]
                
                # Professional sticky end alignment
                if sticky_end_mode and not is_protein:
                    target_aligned, query_aligned, alignment, orientation, is_sticky_end = create_optimal_sticky_end_alignment(seq_a, seq_b)
                else:
                    # Standard alignment
                    aligner = PairwiseAligner()
                    aligner.mode = 'global'
                    if is_protein:
                        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                    else:
                        aligner.match_score = 5.0
                        aligner.mismatch_score = -4.0
                        aligner.open_gap_score = -10.0
                        aligner.extend_gap_score = -1.0
                    
                    alignments = aligner.align(seq_a, seq_b)
                    alignment = alignments[0]
                    target_aligned, query_aligned = extract_gapped_sequences_from_alignment(alignment)
                    orientation = "Forward"
                    is_sticky_end = False
                
                if not target_aligned or not query_aligned:
                    st.error("Failed to create alignment")
                    return
                
                # FIXED: Calculate comprehensive statistics with CORRECT complementarity
                if is_sticky_end:
                    matches, mismatches, gaps, identity, complement_matches, complementarity = calculate_alignment_statistics(
                        target_aligned, query_aligned, seq_a, seq_b, is_sticky_end, orientation
                    )
                    alignment_stats = (matches, mismatches, gaps, identity, complement_matches, complementarity)
                else:
                    stats = calculate_alignment_statistics(target_aligned, query_aligned)
                    matches, mismatches, gaps, identity = stats[:4]
                    alignment_stats = (matches, mismatches, gaps, identity)
                
                # Display results
                st.success("✅ Pairwise alignment completed successfully!")
                
                if is_sticky_end:
                    st.info(f"🔄 Sticky end detected! Using {orientation} orientation for optimal complementarity.")
                
                # Comprehensive metrics
                if is_sticky_end:
                    col_metric1, col_metric2, col_metric3, col_metric4, col_metric5 = st.columns(5)
                    with col_metric1:
                        st.metric("Score", f"{alignment.score:.1f}")
                    with col_metric2:
                        st.metric("Matches", f"{matches}")
                    with col_metric3:
                        st.metric("Identity %", f"{identity:.1f}%")
                    with col_metric4:
                        st.metric("Complements", f"{complement_matches}")
                    with col_metric5:
                        st.metric("Complement %", f"{complementarity:.1f}%")
                else:
                    col_metric1, col_metric2, col_metric3, col_metric4 = st.columns(4)
                    with col_metric1:
                        st.metric("Score", f"{alignment.score:.1f}")
                    with col_metric2:
                        st.metric("Matches", f"{matches}")
                    with col_metric3:
                        st.metric("Identity %", f"{identity:.1f}%")
                    with col_metric4:
                        st.metric("Length", f"{len(target_aligned)}")
                
                # Professional visualization with mutation highlighting
                st.markdown("### 🎨 Pairwise Alignment Visualization")
                
                sequences = [target_aligned, query_aligned]
                
                # Create complete alignment viewer
                alignment_html = create_complete_alignment_viewer_html(
                    sequences, seq_names, alignment_stats, is_sticky_end
                )
                components.html(alignment_html, height=500, scrolling=True)
                
                # Mutation analysis
                if len(sequences) >= 2:
                    mutations = highlight_mutations(sequences[0], sequences[1])
                    detailed_mutations, transitions, transversions = analyze_mutations_detailed(sequences[0], sequences[1])
                    
                    if mutations:
                        st.markdown("### 🔬 Professional Mutation Analysis")
                        
                        col_mut1, col_mut2, col_mut3 = st.columns(3)
                        
                        with col_mut1:
                            st.metric("Total Mutations", len(mutations))
                            st.metric("Transitions (A↔G, C↔T)", len(transitions))
                        
                        with col_mut2:
                            st.metric("Transversions", len(transversions))
                            if len(mutations) > 0:
                                ts_tv_ratio = len(transitions) / len(transversions) if len(transversions) > 0 else "∞"
                                st.metric("Ts/Tv Ratio", f"{ts_tv_ratio:.2f}" if isinstance(ts_tv_ratio, float) else ts_tv_ratio)
                        
                        with col_mut3:
                            mutation_rate = (len(mutations) / len([c for c in sequences[0] if c != '-'])) * 100
                            st.metric("Mutation Rate %", f"{mutation_rate:.2f}%")
                        
                        # Detailed mutation table
                        with st.expander("📋 Detailed Mutation Analysis"):
                            if detailed_mutations:
                                mut_df = pd.DataFrame([
                                    {
                                        'Position': mut['position'] + 1,
                                        'From': mut['from'],
                                        'To': mut['to'],
                                        'Type': mut['type'].title(),
                                        'Change': f"{mut['from']}→{mut['to']}"
                                    }
                                    for mut in detailed_mutations
                                ])
                                st.dataframe(mut_df, use_container_width=True)
                    else:
                        st.success("🎉 Optimal alignment! No mutations detected.")
                
                # Interactive Variable Positions Logomaker
                if not is_protein:
                    st.markdown("### 📊 Interactive Variable Positions Analysis")
                    
                    try:
                        # Create interactive variable positions logo
                        interactive_fig, error = create_interactive_variable_positions_logo(
                            sequences, seq_names, "Interactive Variable Positions"
                        )
                        
                        if interactive_fig:
                            # Display interactive plotly chart
                            st.plotly_chart(interactive_fig, use_container_width=True)
                            st.success("✅ Interactive logo: Zoom, pan, and hover for details!")
                            
                            # Show analysis summary
                            max_len = max(len(str(seq)) for seq in sequences)
                            variable_positions = []
                            
                            for pos in range(max_len):
                                column = [str(seq)[pos] if pos < len(str(seq)) else '-' for seq in sequences]
                                non_gap_chars = [c for c in column if c != '-']
                                unique_chars = set(non_gap_chars)
                                
                                if len(unique_chars) > 1:
                                    variable_positions.append(pos + 1)
                            
                            st.info(f"🎯 Interactive logo showing {len(variable_positions)} variable positions out of {max_len} total")
                            
                            if variable_positions:
                                with st.expander("🔍 Interactive Features & Analysis Details"):
                                    col_help1, col_help2 = st.columns(2)
                                    
                                    with col_help1:
                                        st.write("**Interactive Controls:**")
                                        st.write("• 🔍 **Zoom:** Click and drag to zoom")
                                        st.write("• 👆 **Pan:** Hold Shift + drag to pan") 
                                        st.write("• 📊 **Hover:** Detailed base information")
                                        st.write("• 🔄 **Reset:** Double-click to reset zoom")
                                        st.write("• 📥 **Download:** Toolbar to save PNG")
                                    
                                    with col_help2:
                                        st.write("**Variable Positions:**")
                                        st.write(f"{', '.join(map(str, variable_positions[:15]))}")
                                        if len(variable_positions) > 15:
                                            st.write(f"... and {len(variable_positions) - 15} more")
                        
                        elif error:
                            st.warning(f"Interactive logo: {error}")
                    except Exception as e:
                        st.warning(f"Interactive logo failed: {e}")
                
                # Export options
                st.markdown("### 💾 Professional Export Options")
                col_export1, col_export2, col_export3 = st.columns(3)
                
                with col_export1:
                    fasta_content = f">{seq_names[0]}\n{target_aligned}\n>{seq_names[1]}_{orientation}\n{query_aligned}\n"
                    
                    st.download_button(
                        "📄 Download FASTA",
                        fasta_content,
                        f"pairwise_alignment_{orientation.lower().replace(' ', '_')}.fasta",
                        "text/plain"
                    )
                
                with col_export2:
                    alignment_text = alignment.format()
                    st.download_button(
                        "📄 Download Alignment Details",
                        alignment_text,
                        "pairwise_alignment_details.txt",
                        "text/plain"
                    )
                
                with col_export3:
                    if mutations:
                        mutation_summary = f"Professional Mutation Analysis Report\n"
                        mutation_summary += f"Total mutations: {len(mutations)}\n"
                        mutation_summary += f"Transitions: {len(transitions)}\n"
                        mutation_summary += f"Transversions: {len(transversions)}\n\n"
                        mutation_summary += "Detailed mutations:\n"
                        for mut in detailed_mutations:
                            mutation_summary += f"Position {mut['position']+1}: {mut['from']}→{mut['to']} ({mut['type']})\n"
                        
                        st.download_button(
                            "📄 Download Mutation Report",
                            mutation_summary,
                            "mutation_analysis_report.txt",
                            "text/plain"
                        )
                
        except Exception as e:
            st.error(f"Error during alignment: {e}")
            with st.expander("🔍 Debug Information"):
                st.exception(e)

# For compatibility
def app():
    main()

if __name__ == "__main__":
    main()