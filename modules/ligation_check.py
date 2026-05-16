# modules/ligation_check.py
"""
G-Synth Compatible Synthetic Insert Ligation Check Module
Final corrected version with all fixes applied
- Fixed text area height error (minimum 68px requirement)
- Manual name input (optional)
- Improved layout with visualization in results section
- Enhanced SnapGene-style visualizations
- Complete G-Synth compatibility
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import io
import re
import json
import base64
import math
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any

# Bio_utils imports with comprehensive error handling
try:
    from utils.bio_utils import (
        clean_dna_sequence,
        reverse_complement,
        calculate_tm_consensus,
        validate_dna_sequence,
        calculate_gc
    )
    BIO_UTILS_AVAILABLE = True
except ImportError as e:
    st.error(f"❌ Bio utils import failed: {e}")
    BIO_UTILS_AVAILABLE = False

# ═══════════════════════════════════════════════════════════════════════════════
# G-SYNTH ENZYME DEFINITIONS AND DATABASE
# ═══════════════════════════════════════════════════════════════════════════════

# G-Synth compatible enzyme pairs (authoritative for ligation)
GSYNTH_ENZYME_PAIRS = {
    "NdeI / XhoI": {"forward_overhang": "TA", "reverse_overhang": "TCGA"},
    "NdeI / EcoRI": {"forward_overhang": "TA", "reverse_overhang": "AATT"},
    "NdeI / BamHI": {"forward_overhang": "TA", "reverse_overhang": "GATC"},
    "BamHI / EcoRI": {"forward_overhang": "GATC", "reverse_overhang": "AATT"},
    "BamHI / XhoI": {"forward_overhang": "GATC", "reverse_overhang": "TCGA"},
    "BamHI / HindIII": {"forward_overhang": "GATC", "reverse_overhang": "AGCT"},
    "EcoRI / HindIII": {"forward_overhang": "AATT", "reverse_overhang": "AGCT"},
    "EcoRI / XbaI": {"forward_overhang": "AATT", "reverse_overhang": "CTAG"},
    "SalI / XbaI": {"forward_overhang": "TCGAC", "reverse_overhang": "TCTAG"},
    "SalI / NotI": {"forward_overhang": "TCGAC", "reverse_overhang": "GGCC"},
    "XhoI / NotI": {"forward_overhang": "TCGA", "reverse_overhang": "GGCC"},
    "KpnI / XbaI": {"forward_overhang": "GGTAC", "reverse_overhang": "CTAG"},
    "SpeI / XbaI": {"forward_overhang": "CTAG", "reverse_overhang": "CTAG"},
    "NcoI / XhoI": {"forward_overhang": "CATG", "reverse_overhang": "TCGA"}
}

# Enhanced enzyme database with comprehensive information
GSYNTH_ENZYME_DATABASE = {
    "NdeI": {
        "recognition": "CATATG", "cut_position": 2, "overhang": "TATG", "overhang_for_ligation": "TA",
        "type": "Type II", "source": "Neisseria denitrificans", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": False,
        "description": "Expression cloning enzyme", "color": "#FF6B35",
        "applications": ["Protein expression", "N-terminal cloning", "His-tag vectors"]
    },
    "XhoI": {
        "recognition": "CTCGAG", "cut_position": 1, "overhang": "TCGA", "overhang_for_ligation": "TCGA",
        "type": "Type II", "source": "Xanthomonas holcicola", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": False,
        "description": "Popular cloning enzyme", "color": "#4ECDC4",
        "applications": ["General cloning", "C-terminal cloning", "Gateway cloning"]
    },
    "EcoRI": {
        "recognition": "GAATTC", "cut_position": 1, "overhang": "AATT", "overhang_for_ligation": "AATT",
        "type": "Type II", "source": "Escherichia coli", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": True,
        "description": "Classic restriction enzyme", "color": "#45B7D1",
        "applications": ["General cloning", "Library construction", "Mapping"]
    },
    "BamHI": {
        "recognition": "GGATCC", "cut_position": 1, "overhang": "GATC", "overhang_for_ligation": "GATC",
        "type": "Type II", "source": "Bacillus amyloliquefaciens", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": True,
        "description": "Very popular cloning enzyme", "color": "#96CEB4",
        "applications": ["General cloning", "Subcloning", "Vector construction"]
    },
    "HindIII": {
        "recognition": "AAGCTT", "cut_position": 1, "overhang": "AGCT", "overhang_for_ligation": "AGCT",
        "type": "Type II", "source": "Haemophilus influenzae", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": True,
        "description": "Traditional restriction enzyme", "color": "#FFEAA7",
        "applications": ["General cloning", "Mapping", "Southern blot"]
    },
    "SalI": {
        "recognition": "GTCGAC", "cut_position": 1, "overhang": "TCGA", "overhang_for_ligation": "TCGAC",
        "type": "Type II", "source": "Streptomyces albus", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": False,
        "description": "Compatible with XhoI", "color": "#DDA0DD",
        "applications": ["Directional cloning", "Compatible cloning", "Subcloning"]
    },
    "XbaI": {
        "recognition": "TCTAGA", "cut_position": 1, "overhang": "CTAG", "overhang_for_ligation": "TCTAG",
        "type": "Type II", "source": "Xanthomonas badrii", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": True,
        "description": "BioBrick standard enzyme", "color": "#FF7675",
        "applications": ["BioBrick assembly", "Synthetic biology", "Modular cloning"]
    },
    "NotI": {
        "recognition": "GCGGCCGC", "cut_position": 2, "overhang": "GGCC", "overhang_for_ligation": "GGCC",
        "type": "Type II", "source": "Nocardia otitidis", "temperature": 37,
        "buffer": "NEB Buffer 3.1", "heat_inactivation": False, "star_activity": False,
        "description": "Rare 8-base cutter", "color": "#A29BFE",
        "applications": ["Large fragment cloning", "BAC cloning", "Genomic libraries"]
    },
    "NcoI": {
        "recognition": "CCATGG", "cut_position": 1, "overhang": "CATG", "overhang_for_ligation": "CATG",
        "type": "Type II", "source": "Nocardia corallina", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": False,
        "description": "Alternative to NdeI", "color": "#FD79A8",
        "applications": ["Protein expression", "Alternative start site", "Expression vectors"]
    },
    "KpnI": {
        "recognition": "GGTACC", "cut_position": 5, "overhang": "GGTAC", "overhang_for_ligation": "GGTAC",
        "type": "Type II", "source": "Klebsiella pneumoniae", "temperature": 37,
        "buffer": "NEB Buffer 1.1", "heat_inactivation": True, "star_activity": False,
        "description": "Creates 3' overhang", "color": "#FDCB6E",
        "applications": ["3' overhang cloning", "Specific orientation", "Vector linearization"]
    },
    "SpeI": {
        "recognition": "ACTAGT", "cut_position": 1, "overhang": "CTAG", "overhang_for_ligation": "CTAG",
        "type": "Type II", "source": "Sphaerotilus natans", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": False,
        "description": "Compatible with XbaI", "color": "#E17055",
        "applications": ["Compatible cloning", "BioBrick assembly", "Fusion cloning"]
    },
    "BglII": {
        "recognition": "AGATCT", "cut_position": 1, "overhang": "GATC", "overhang_for_ligation": "GATC",
        "type": "Type II", "source": "Bacillus globigii", "temperature": 37,
        "buffer": "NEB Buffer 2.1", "heat_inactivation": True, "star_activity": False,
        "description": "Compatible with BamHI", "color": "#00B894",
        "applications": ["Compatible cloning", "Subcloning", "Fusion proteins"]
    }
}

# ═══════════════════════════════════════════════════════════════════════════════
# G-SYNTH ANALYSIS FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def gsynth_find_restriction_sites(sequence: str, enzyme_name: str) -> List[Dict]:
    """Find all restriction sites for a given enzyme in a sequence."""
    if enzyme_name not in GSYNTH_ENZYME_DATABASE:
        return []
    
    enzyme_data = GSYNTH_ENZYME_DATABASE[enzyme_name]
    recognition = enzyme_data["recognition"]
    sites = []
    
    start = 0
    while True:
        pos = sequence.find(recognition, start)
        if pos == -1:
            break
        
        sites.append({
            "enzyme": enzyme_name,
            "position": pos,
            "recognition": recognition,
            "sequence": sequence[pos:pos+len(recognition)],
            "cut_position": pos + enzyme_data["cut_position"],
            "overhang": enzyme_data["overhang"],
            "context_5": sequence[max(0, pos-10):pos],
            "context_3": sequence[pos+len(recognition):pos+len(recognition)+10]
        })
        start = pos + 1
    
    return sites

def gsynth_simulate_digestion_enhanced(vector_seq: str, enzyme_pair: str) -> Dict:
    """Enhanced digestion simulation with detailed analysis."""
    try:
        enzymes = enzyme_pair.split(" / ")
        if len(enzymes) != 2:
            return {"error": "Invalid enzyme pair format"}
        
        forward_enzyme, reverse_enzyme = enzymes
        
        # Find all restriction sites
        forward_sites = gsynth_find_restriction_sites(vector_seq, forward_enzyme)
        reverse_sites = gsynth_find_restriction_sites(vector_seq, reverse_enzyme)
        
        if not forward_sites:
            return {"error": f"No {forward_enzyme} sites found in vector"}
        if not reverse_sites:
            return {"error": f"No {reverse_enzyme} sites found in vector"}
        
        # Use the first occurrence of each site
        forward_site = forward_sites[0]
        reverse_site = reverse_sites[0]
        
        forward_cut = forward_site["cut_position"]
        reverse_cut = reverse_site["cut_position"]
        
        # Determine cut order and extract fragments
        if forward_cut < reverse_cut:
            # Forward site comes first
            linearized_vector = vector_seq[:forward_cut] + vector_seq[reverse_cut:]
            excised_fragment = vector_seq[forward_cut:reverse_cut]
        else:
            # Reverse site comes first (circular)
            linearized_vector = vector_seq[reverse_cut:] + vector_seq[:forward_cut]
            excised_fragment = vector_seq[forward_cut:] + vector_seq[:reverse_cut]
        
        return {
            "success": True,
            "original_vector": vector_seq,
            "linearized_vector": linearized_vector,
            "excised_fragment": excised_fragment,
            "forward_sites": forward_sites,
            "reverse_sites": reverse_sites,
            "cut_positions": [forward_cut, reverse_cut],
            "enzyme_pair": enzyme_pair,
            "fragments": {
                "vector_length": len(linearized_vector),
                "fragment_length": len(excised_fragment),
                "total_recovered": len(linearized_vector) + len(excised_fragment)
            }
        }
        
    except Exception as e:
        return {"error": f"Digestion simulation failed: {str(e)}"}

def gsynth_validate_sticky_ends_enhanced(forward_seq: str, reverse_seq: str, enzyme_pair: str) -> Dict:
    """Enhanced sticky end validation with detailed analysis."""
    try:
        pair_info = GSYNTH_ENZYME_PAIRS.get(enzyme_pair, {})
        if not pair_info:
            return {"error": f"Unknown enzyme pair: {enzyme_pair}"}
        
        forward_overhang = pair_info.get("forward_overhang", "")
        reverse_overhang = pair_info.get("reverse_overhang", "")
        
        # Check forward end
        forward_start = forward_seq[:len(forward_overhang)] if len(forward_seq) >= len(forward_overhang) else forward_seq
        forward_valid = forward_start == forward_overhang
        forward_similarity = (sum(a == b for a, b in zip(forward_start, forward_overhang)) / len(forward_overhang) * 100) if forward_overhang else 100
        
        # Check reverse end
        reverse_start = reverse_seq[:len(reverse_overhang)] if len(reverse_seq) >= len(reverse_overhang) else reverse_seq
        reverse_valid = reverse_start == reverse_overhang
        reverse_similarity = (sum(a == b for a, b in zip(reverse_start, reverse_overhang)) / len(reverse_overhang) * 100) if reverse_overhang else 100
        
        return {
            "overall_valid": forward_valid and reverse_valid,
            "forward_analysis": {
                "valid": forward_valid,
                "expected": forward_overhang,
                "found": forward_start,
                "similarity": forward_similarity,
                "enzyme": enzyme_pair.split(" / ")[0]
            },
            "reverse_analysis": {
                "valid": reverse_valid,
                "expected": reverse_overhang,
                "found": reverse_start,
                "similarity": reverse_similarity,
                "enzyme": enzyme_pair.split(" / ")[1]
            }
        }
        
    except Exception as e:
        return {"error": f"Sticky end validation failed: {str(e)}"}

def gsynth_check_hybridization_enhanced(forward_seq: str, reverse_seq: str, enzyme_pair: str) -> Dict:
    """Enhanced hybridization check with sequence alignment."""
    try:
        pair_info = GSYNTH_ENZYME_PAIRS.get(enzyme_pair, {})
        if not pair_info:
            return {"error": f"Unknown enzyme pair: {enzyme_pair}"}
        
        forward_overhang = pair_info.get("forward_overhang", "")
        reverse_overhang = pair_info.get("reverse_overhang", "")
        
        # Extract internal sequences (without overhangs)
        forward_internal = forward_seq[len(forward_overhang):] if len(forward_seq) > len(forward_overhang) else ""
        reverse_internal = reverse_seq[len(reverse_overhang):] if len(reverse_seq) > len(reverse_overhang) else ""
        
        # Get reverse complement of reverse internal sequence
        reverse_complement_seq = reverse_complement(reverse_internal) if BIO_UTILS_AVAILABLE else reverse_internal[::-1]
        
        # Check if they match (perfect hybridization)
        min_length = min(len(forward_internal), len(reverse_complement_seq))
        if min_length == 0:
            match_percentage = 0
            perfect_match = False
        else:
            matches = sum(a == b for a, b in zip(forward_internal[:min_length], reverse_complement_seq[:min_length]))
            match_percentage = (matches / min_length) * 100
            perfect_match = forward_internal == reverse_complement_seq
        
        return {
            "hybridization_success": perfect_match,
            "match_percentage": match_percentage,
            "forward_internal": forward_internal,
            "reverse_internal": reverse_internal,
            "reverse_complement": reverse_complement_seq,
            "alignment_length": min_length,
            "matches": matches if min_length > 0 else 0
        }
        
    except Exception as e:
        return {"error": f"Hybridization check failed: {str(e)}"}

# ═══════════════════════════════════════════════════════════════════════════════
# ENHANCED SNAPGENE-STYLE VISUALIZATION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def create_enhanced_plasmid_map(vector_seq: str, cut_sites: List[int], enzyme_pair: str, 
                               insert_seq: str = None, success: bool = False,
                               plasmid_name: str = "Vector", insert_name: str = "Insert") -> go.Figure:
    """Create enhanced SnapGene-style plasmid map with accurate positioning."""
    
    fig = go.Figure()
    
    # Calculate plasmid properties
    if success and insert_seq:
        total_length = len(vector_seq) + len(insert_seq)
        display_name = f"{plasmid_name} + {insert_name}"
    else:
        total_length = len(vector_seq)
        display_name = plasmid_name
    
    # Enhanced circular coordinates
    theta = np.linspace(0, 2*np.pi, 300)
    outer_radius = 1.0
    inner_radius = 0.65
    
    if success and insert_seq:
        # Successful ligation visualization
        vector_proportion = len(vector_seq) / total_length
        
        # Vector backbone arc
        theta_vector = theta[:int(vector_proportion * 300)]
        x_vector_outer = outer_radius * np.cos(theta_vector)
        y_vector_outer = outer_radius * np.sin(theta_vector)
        x_vector_inner = inner_radius * np.cos(theta_vector)
        y_vector_inner = inner_radius * np.sin(theta_vector)
        
        # Insert arc
        theta_insert = theta[int(vector_proportion * 300):]
        x_insert_outer = outer_radius * np.cos(theta_insert)
        y_insert_outer = outer_radius * np.sin(theta_insert)
        x_insert_inner = inner_radius * np.cos(theta_insert)
        y_insert_inner = inner_radius * np.sin(theta_insert)
        
        # Vector backbone with enhanced styling
        fig.add_trace(go.Scatter(
            x=np.concatenate([x_vector_outer, x_vector_inner[::-1], [x_vector_outer[0]]]),
            y=np.concatenate([y_vector_outer, y_vector_inner[::-1], [y_vector_outer[0]]]),
            fill='toself',
            fillcolor='rgba(70, 130, 180, 0.85)',
            line=dict(color='#2C5282', width=4),
            name=f'{plasmid_name} ({len(vector_seq):,} bp)',
            hovertemplate=f'{plasmid_name}<br>Length: {len(vector_seq):,} bp<br>GC: {calculate_gc(vector_seq) if BIO_UTILS_AVAILABLE else "N/A":.1f}%<extra></extra>',
            showlegend=True
        ))
        
        # Insert with contrasting color
        fig.add_trace(go.Scatter(
            x=np.concatenate([x_insert_outer, x_insert_inner[::-1], [x_insert_outer[0]]]),
            y=np.concatenate([y_insert_outer, y_insert_inner[::-1], [y_insert_outer[0]]]),
            fill='toself',
            fillcolor='rgba(34, 197, 94, 0.85)',
            line=dict(color='#16A34A', width=4),
            name=f'{insert_name} ({len(insert_seq):,} bp)',
            hovertemplate=f'{insert_name}<br>Length: {len(insert_seq):,} bp<br>GC: {calculate_gc(insert_seq) if BIO_UTILS_AVAILABLE else "N/A":.1f}%<extra></extra>',
            showlegend=True
        ))
        
        # Enhanced ligation junctions
        junction_angles = [theta_vector[-1] if len(theta_vector) > 0 else 0, 
                          theta_insert[0] if len(theta_insert) > 0 else np.pi]
        for i, angle in enumerate(junction_angles):
            x_junction = (outer_radius + 0.08) * np.cos(angle)
            y_junction = (outer_radius + 0.08) * np.sin(angle)
            
            fig.add_trace(go.Scatter(
                x=[x_junction], y=[y_junction],
                mode='markers',
                marker=dict(
                    size=22,
                    color='#F59E0B',
                    symbol='star',
                    line=dict(color='white', width=3)
                ),
                name='Ligation Site' if i == 0 else None,
                showlegend=i == 0,
                hovertemplate='Ligation Junction<br>T4 DNA Ligase binding<extra></extra>'
            ))
        
    else:
        # Original circular vector
        x_outer = outer_radius * np.cos(theta)
        y_outer = outer_radius * np.sin(theta)
        x_inner = inner_radius * np.cos(theta)
        y_inner = inner_radius * np.sin(theta)
        
        fig.add_trace(go.Scatter(
            x=np.concatenate([x_outer, x_inner[::-1], [x_outer[0]]]),
            y=np.concatenate([y_outer, y_inner[::-1], [y_outer[0]]]),
            fill='toself',
            fillcolor='rgba(70, 130, 180, 0.85)',
            line=dict(color='#2C5282', width=4),
            name=f'{plasmid_name} ({len(vector_seq):,} bp)',
            hovertemplate=f'{plasmid_name}<br>Length: {len(vector_seq):,} bp<br>GC: {calculate_gc(vector_seq) if BIO_UTILS_AVAILABLE else "N/A":.1f}%<extra></extra>',
            showlegend=True
        ))
    
    # Enhanced restriction enzyme sites with improved positioning
    if cut_sites and enzyme_pair:
        enzymes = enzyme_pair.split(" / ")
        enzyme_colors = ['#EF4444', '#F59E0B']
        
        for i, (site_pos, enzyme) in enumerate(zip(cut_sites, enzymes)):
            angle = (site_pos / len(vector_seq)) * 2 * np.pi
            
            # Enzyme marker with improved positioning
            marker_radius = outer_radius + 0.20
            x_site = marker_radius * np.cos(angle)
            y_site = marker_radius * np.sin(angle)
            
            fig.add_trace(go.Scatter(
                x=[x_site], y=[y_site],
                mode='markers+text',
                marker=dict(
                    size=28,
                    color=enzyme_colors[i % len(enzyme_colors)],
                    symbol='triangle-up',
                    line=dict(color='white', width=3)
                ),
                text=[enzyme],
                textposition='middle center',
                textfont=dict(size=11, color='white', family='Arial Black'),
                name=f'{enzyme} Site',
                hovertemplate=f'{enzyme} Recognition Site<br>Position: {site_pos:,} bp<br>Recognition: {GSYNTH_ENZYME_DATABASE.get(enzyme, {}).get("recognition", "Unknown")}<br>Overhang: {GSYNTH_ENZYME_DATABASE.get(enzyme, {}).get("overhang", "Unknown")}<extra></extra>'
            ))
            
            # Enhanced cutting arrow with precise positioning
            arrow_start_radius = outer_radius + 0.15
            arrow_end_radius = outer_radius + 0.02
            x_arrow_start = arrow_start_radius * np.cos(angle)
            y_arrow_start = arrow_start_radius * np.sin(angle)
            x_arrow_end = arrow_end_radius * np.cos(angle)
            y_arrow_end = arrow_end_radius * np.sin(angle)
            
            fig.add_annotation(
                x=x_arrow_end, y=y_arrow_end,
                ax=x_arrow_start, ay=y_arrow_start,
                arrowhead=3, arrowsize=1.8, arrowwidth=4,
                arrowcolor=enzyme_colors[i % len(enzyme_colors)],
                showarrow=True
            )
            
            # Add enzyme name label with improved positioning
            label_radius = outer_radius + 0.35
            x_label = label_radius * np.cos(angle)
            y_label = label_radius * np.sin(angle)
            
            # Calculate text angle for better readability
            text_angle = angle * 180 / np.pi
            if text_angle > 90 and text_angle < 270:
                text_angle += 180  # Flip text to avoid upside-down reading
            
            fig.add_annotation(
                x=x_label, y=y_label,
                text=f"{enzyme}<br>{site_pos:,} bp",
                showarrow=False,
                font=dict(size=10, color=enzyme_colors[i % len(enzyme_colors)], family='Arial Black'),
                align="center",
                textangle=text_angle if abs(text_angle) < 45 else 0,
                bgcolor="rgba(255, 255, 255, 0.9)",
                bordercolor=enzyme_colors[i % len(enzyme_colors)],
                borderwidth=1,
                borderpad=3
            )
    
    # Enhanced scale markers with better positioning
    scale_positions = [0, np.pi/2, np.pi, 3*np.pi/2]
    scale_labels = ['0', f'{total_length//4:,}', f'{total_length//2:,}', f'{3*total_length//4:,}']
    
    for pos, label in zip(scale_positions, scale_labels):
        x_scale = (inner_radius - 0.15) * np.cos(pos)
        y_scale = (inner_radius - 0.15) * np.sin(pos)
        
        fig.add_annotation(
            x=x_scale, y=y_scale,
            text=label,
            showarrow=False,
            font=dict(size=10, color='#374151', family='Arial'),
            align="center"
        )
        
        # Enhanced tick marks
        x_tick_inner = (inner_radius - 0.08) * np.cos(pos)
        y_tick_inner = (inner_radius - 0.08) * np.sin(pos)
        x_tick_outer = (inner_radius + 0.08) * np.cos(pos)
        y_tick_outer = (inner_radius + 0.08) * np.sin(pos)
        
        fig.add_shape(
            type="line",
            x0=x_tick_inner, y0=y_tick_inner,
            x1=x_tick_outer, y1=y_tick_outer,
            line=dict(color='#6B7280', width=3)
        )
    
    # Enhanced central information with custom names and insert size
    center_text = f"<b>{display_name}</b><br>{total_length:,} bp"
    if success and insert_seq:
        center_text += f"<br>Insert: {len(insert_seq):,} bp"
    if enzyme_pair:
        center_text += f"<br><i>{enzyme_pair}</i>"
    
    fig.add_annotation(
        x=0, y=0,
        text=center_text,
        showarrow=False,
        font=dict(size=14, color='#1F2937', family='Arial'),
        align="center",
        bgcolor="rgba(255, 255, 255, 0.95)",
        bordercolor="#D1D5DB",
        borderwidth=2,
        borderpad=12
    )
    
    # Enhanced title with success indication
    title_text = f"✅ {display_name} - Recombinant Plasmid" if success else f"🧬 {plasmid_name} - Restriction Map"
    title_color = '#10B981' if success else '#3B82F6'
    
    fig.update_layout(
        title=dict(
            text=title_text,
            font=dict(size=20, color=title_color, family='Arial Black'),
            x=0.5, y=0.95
        ),
        xaxis=dict(visible=False, range=[-1.6, 1.6]),
        yaxis=dict(visible=False, range=[-1.6, 1.6]),
        plot_bgcolor='rgba(248, 250, 252, 0.98)',
        paper_bgcolor='white',
        height=650,
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.98,
            xanchor="left",
            x=1.02,
            bgcolor="rgba(255, 255, 255, 0.95)",
            bordercolor="#D1D5DB",
            borderwidth=1
        )
    )
    
    # Ensure equal aspect ratio
    fig.update_xaxes(scaleanchor="y", scaleratio=1)
    
    return fig

def create_enhanced_junction_analysis(forward_seq: str, reverse_seq: str, enzyme_pair: str,
                                     validation_result: Dict, insert_name: str = "Insert") -> go.Figure:
    """Create enhanced junction analysis with molecular detail - FIXED VERSION."""
    
    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=[
            f'5\' → 3\' Forward Strand ({insert_name})',
            f'5\' → 3\' Reverse Strand ({insert_name})', 
            'Hybridization Analysis'
        ],
        vertical_spacing=0.15,
        row_heights=[0.3, 0.3, 0.4]
    )
    
    # Get enzyme information
    pair_info = GSYNTH_ENZYME_PAIRS.get(enzyme_pair, {})
    forward_overhang = pair_info.get("forward_overhang", "")
    reverse_overhang = pair_info.get("reverse_overhang", "")
    
    # Enhanced color scheme for bases (SnapGene-style)
    base_colors = {
        'A': '#FF6B6B',  # Red
        'T': '#4ECDC4',  # Teal  
        'G': '#45B7D1',  # Blue
        'C': '#96CEB4',  # Green
        'N': '#95A5A6'   # Gray
    }
    
    # Forward strand visualization
    for i, base in enumerate(forward_seq):
        is_overhang = i < len(forward_overhang)
        base_color = base_colors.get(base, '#95A5A6')
        marker_size = 32 if is_overhang else 28
        opacity = 1.0 if is_overhang else 0.85
        
        fig.add_trace(go.Scatter(
            x=[i], y=[1],
            mode='markers+text',
            marker=dict(
                size=marker_size,
                color=base_color,
                opacity=opacity,
                line=dict(color='white', width=2)
            ),
            text=[base],
            textposition='middle center',
            textfont=dict(size=12, color='white', family='Arial Black'),
            name='Forward' if i == 0 else None,
            showlegend=i == 0,
            hovertemplate=f'Position: {i}<br>Base: {base}<br>Type: {"Sticky End" if is_overhang else "Internal"}<extra></extra>'
        ), row=1, col=1)
    
    # Reverse strand visualization
    for i, base in enumerate(reverse_seq):
        is_overhang = i < len(reverse_overhang)
        base_color = base_colors.get(base, '#95A5A6')
        marker_size = 32 if is_overhang else 28
        opacity = 1.0 if is_overhang else 0.85
        
        fig.add_trace(go.Scatter(
            x=[i], y=[1],
            mode='markers+text',
            marker=dict(
                size=marker_size,
                color=base_color,
                opacity=opacity,
                line=dict(color='white', width=2)
            ),
            text=[base],
            textposition='middle center',
            textfont=dict(size=12, color='white', family='Arial Black'),
            name='Reverse' if i == 0 else None,
            showlegend=i == 0,
            hovertemplate=f'Position: {i}<br>Base: {base}<br>Type: {"Sticky End" if is_overhang else "Internal"}<extra></extra>'
        ), row=2, col=1)
    
    # Enhanced overhang highlighting
    if len(forward_overhang) > 0:
        fig.add_shape(
            type="rect",
            x0=-0.5, y0=0.6, x1=len(forward_overhang)-0.5, y1=1.4,
            line=dict(color="#9C27B0", width=3),
            fillcolor="rgba(156, 39, 176, 0.25)",
            row=1, col=1
        )
        
        fig.add_annotation(
            x=len(forward_overhang)/2 - 0.5, y=1.7,
            text=f"5' Overhang: {forward_overhang}",
            showarrow=False,
            font=dict(size=12, color="#9C27B0", family='Arial Black'),
            row=1, col=1
        )
    
    if len(reverse_overhang) > 0:
        fig.add_shape(
            type="rect",
            x0=-0.5, y0=0.6, x1=len(reverse_overhang)-0.5, y1=1.4,
            line=dict(color="#E91E63", width=3),
            fillcolor="rgba(233, 30, 99, 0.25)",
            row=2, col=1
        )
        
        fig.add_annotation(
            x=len(reverse_overhang)/2 - 0.5, y=1.7,
            text=f"5' Overhang: {reverse_overhang}",
            showarrow=False,
            font=dict(size=12, color="#E91E63", family='Arial Black'),
            row=2, col=1
        )
    
    # Hybridization analysis in third panel
    if "forward_internal" in validation_result and "reverse_complement" in validation_result:
        forward_internal = validation_result["forward_internal"]
        reverse_comp = validation_result["reverse_complement"]
        match_percentage = validation_result.get("match_percentage", 0)
        
        max_len = max(len(forward_internal), len(reverse_comp))
        
        if max_len > 0:
            # Forward strand
            for i, base in enumerate(forward_internal):
                fig.add_trace(go.Scatter(
                    x=[i], y=[2],
                    mode='markers+text',
                    marker=dict(size=28, color=base_colors.get(base, '#95A5A6')),
                    text=[base],
                    textposition='middle center',
                    textfont=dict(size=12, color='white', family='Arial Black'),
                    showlegend=False,
                    hovertemplate=f'Forward Internal: {base}<br>Position: {i}<extra></extra>'
                ), row=3, col=1)
            
            # Reverse complement strand
            for i, base in enumerate(reverse_comp):
                matches = i < len(forward_internal) and base == forward_internal[i]
                border_color = '#10B981' if matches else '#EF4444'
                
                fig.add_trace(go.Scatter(
                    x=[i], y=[1],
                    mode='markers+text',
                    marker=dict(
                        size=28, 
                        color=base_colors.get(base, '#95A5A6'),
                        line=dict(color=border_color, width=3)
                    ),
                    text=[base],
                    textposition='middle center',
                    textfont=dict(size=12, color='white', family='Arial Black'),
                    showlegend=False,
                    hovertemplate=f'Reverse Complement: {base}<br>Position: {i}<br>Match: {"Yes" if matches else "No"}<extra></extra>'
                ), row=3, col=1)
                
                # Add hydrogen bond lines for matches
                if matches:
                    fig.add_trace(go.Scatter(
                        x=[i, i], y=[1.3, 1.7],
                        mode='lines',
                        line=dict(color='#10B981', width=4, dash='dot'),
                        showlegend=False,
                        hoverinfo='skip'
                    ), row=3, col=1)
        
        # FIXED: Add match percentage annotation using valid Plotly properties
        fig.add_annotation(
            x=max_len/2 if max_len > 0 else 0, y=0.5,
            text=f"Hybridization Match: {match_percentage:.1f}%",
            showarrow=False,
            font=dict(size=16, color='#1F2937', family='Arial Black'),
            bgcolor="lightblue",      # ✅ Valid Plotly property
            bordercolor="#3B82F6",    # ✅ Valid Plotly property  
            borderwidth=2,            # ✅ Valid Plotly property
            borderpad=8,              # ✅ Valid Plotly property
            opacity=0.9,              # ✅ Valid Plotly property
            row=3, col=1
        )
    
    # Add strand direction labels
    fig.add_annotation(
        x=-2, y=1,
        text="5'",
        showarrow=False,
        font=dict(size=14, color='#374151', family='Arial Black'),
        row=1, col=1
    )
    
    fig.add_annotation(
        x=len(forward_seq), y=1,
        text="3'",
        showarrow=False,
        font=dict(size=14, color='#374151', family='Arial Black'),
        row=1, col=1
    )
    
    fig.update_layout(
        height=650,
        title=dict(
            text=f"🔬 Molecular Junction Analysis - {insert_name}",
            font=dict(size=20, color='#1F2937', family='Arial Black'),
            x=0.5
        ),
        showlegend=True,
        plot_bgcolor='rgba(248, 250, 252, 0.98)',
        paper_bgcolor='white'
    )
    
    # Hide axes for cleaner look
    for i in range(1, 4):
        fig.update_xaxes(visible=False, row=i, col=1)
        fig.update_yaxes(visible=False, range=[0.3, 2.5], row=i, col=1)
    
    return fig

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN APPLICATION WITH IMPROVED LAYOUT AND FIXED TEXT AREA HEIGHTS
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    """Main application with enhanced layout, manual name input, and fixed heights."""
    
    # Clean title without subtitle
    st.markdown("""
    <div style="
        background: linear-gradient(135deg, #1E3A8A 0%, #2563EB 50%, #3B82F6 100%);
        padding: 3rem;
        border-radius: 25px;
        text-align: center;
        margin-bottom: 2.5rem;
        color: white;
        box-shadow: 0 15px 50px rgba(37, 99, 235, 0.5);
    ">
        <h1 style="margin: 0; font-size: 3.5rem; font-weight: 900;">
            🧬 Synthetic Insert Ligation Check
        </h1>
    </div>
    """, unsafe_allow_html=True)
    
    # Initialize session state
    if "ligation_data" not in st.session_state:
        st.session_state["ligation_data"] = {
            "digestion_result": None,
            "validation_result": None,
            "analysis_history": []
        }
    
    # Enhanced sidebar with examples
    with st.sidebar:
        st.markdown("### 📋 Quick Examples")
        
        examples = {
            "NdeI/XhoI Expression": {
                "vector": "ATGCATATGAAAGCTTGATATCCTCGAGCTGCAGAATTCGGATCCAAGCTTGCATGC",
                "forward": "TATGAAAGGATCCAAGCTTCTCGA",
                "reverse": "TCGAAAGCTTGGATCCTTTCAT",
                "enzyme": "NdeI / XhoI",
                "plasmid": "pET-28a(+)",
                "insert": "Target Gene"
            },
            "BamHI/EcoRI Classic": {
                "vector": "ATGCGGATCCAAGCTTGATATCGAATTCCTGCAGAATTCGGATCCAAGCTT",
                "forward": "GATCATGAAAGAATTC",
                "reverse": "GAATTCTTTCAT",
                "enzyme": "BamHI / EcoRI",
                "plasmid": "pUC19",
                "insert": "Cloned Fragment"
            },
            "SalI/XbaI BioBrick": {
                "vector": "ATGCGTCGACAAGCTTGATATCTCTAGACTGCAGAATTCGGATCCAAGCTT",
                "forward": "TCGACATGAAAGGTACCTAG",
                "reverse": "TCTAGGTACCTTTCAT",
                "enzyme": "SalI / XbaI",
                "plasmid": "pSB1C3",
                "insert": "BioBrick Part"
            }
        }
        
        for name, data in examples.items():
            if st.button(f"🧬 {name}", use_container_width=True):
                st.session_state.update({
                    "example_vector": data["vector"],
                    "example_forward": data["forward"],
                    "example_reverse": data["reverse"],
                    "example_enzyme": data["enzyme"],
                    "example_plasmid": data["plasmid"],
                    "example_insert": data["insert"]
                })
                st.success(f"Loaded: {name}")
                st.rerun()
        
        st.markdown("---")
        st.markdown("### 🎨 Display Options")
        show_molecular_detail = st.checkbox("Molecular Detail", value=True)
        real_time_analysis = st.checkbox("Real-time Validation", value=True)
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # SECTION 1: INPUT SECTION WITH MANUAL NAME FIELDS
    # ═══════════════════════════════════════════════════════════════════════════════
    
    st.markdown("""
    <div style="
        background: linear-gradient(135deg, #F0F9FF, #E0F2FE);
        padding: 2rem;
        border-radius: 20px;
        margin-bottom: 2rem;
        border: 3px solid #0EA5E9;
        box-shadow: 0 8px 25px rgba(14, 165, 233, 0.2);
    ">
        <h2 style="color: #0C4A6E; margin: 0 0 1.5rem 0; font-size: 1.8rem;">
            📝 Input Section
        </h2>
    """, unsafe_allow_html=True)
    
    # MANUAL NAME INPUT (Optional) - Added as requested
    st.markdown("#### 🏷️ Manual Names (Optional)")
    col_name1, col_name2 = st.columns(2)
    
    with col_name1:
        plasmid_name = st.text_input(
            "Plasmid Name (Optional)",
            value=st.session_state.get("example_plasmid", ""),
            placeholder="e.g., pET-28a(+)",
            help="Enter custom name for your plasmid vector"
        )
    
    with col_name2:
        insert_name = st.text_input(
            "Insert Name (Optional)",
            value=st.session_state.get("example_insert", ""),
            placeholder="e.g., Target Gene",
            help="Enter custom name for your insert fragment"
        )
    
    # Set default names if empty
    if not plasmid_name:
        plasmid_name = "Vector"
    if not insert_name:
        insert_name = "Insert"
    
    st.markdown("---")
    
    # Enzyme selection
    st.markdown("#### ⚔️ Restriction Enzyme Selection")
    col_enz1, col_enz2 = st.columns([2, 1])
    
    with col_enz1:
        enzyme_pair = st.selectbox(
            "Select Enzyme Pair",
            list(GSYNTH_ENZYME_PAIRS.keys()),
            index=list(GSYNTH_ENZYME_PAIRS.keys()).index(st.session_state.get("example_enzyme", "NdeI / XhoI")) if st.session_state.get("example_enzyme") in GSYNTH_ENZYME_PAIRS else 0,
            help="Choose restriction enzymes for digestion"
        )
    
    with col_enz2:
        if enzyme_pair in GSYNTH_ENZYME_PAIRS:
            pair_info = GSYNTH_ENZYME_PAIRS[enzyme_pair]
            enzymes = enzyme_pair.split(" / ")
            
            st.info(f"""
            **Overhangs:**
            - {enzymes[0]}: `{pair_info['forward_overhang']}`
            - {enzymes[1]}: `{pair_info['reverse_overhang']}`
            """)
    
    st.markdown("---")
    
    # Sequence input layout - FIXED TEXT AREA HEIGHTS
    col_input1, col_input2 = st.columns([1.5, 1])
    
    with col_input1:
        st.markdown("#### Vector Sequence")
        vector_sequence = st.text_area(
            "Vector DNA sequence",
            value=st.session_state.get("example_vector", ""),
            height=120,  # ✅ Already correct - above minimum
            placeholder="Enter circular vector sequence containing restriction sites...",
            help="Paste your vector sequence here"
        )
        
        # Real-time vector analysis
        if vector_sequence and real_time_analysis:
            clean_vector = clean_dna_sequence(vector_sequence) if BIO_UTILS_AVAILABLE else vector_sequence.upper()
            
            col_rt1, col_rt2, col_rt3 = st.columns(3)
            with col_rt1:
                st.metric("Length", f"{len(clean_vector):,} bp")
            with col_rt2:
                gc_content = calculate_gc(clean_vector) if BIO_UTILS_AVAILABLE else 0
                st.metric("GC%", f"{gc_content:.1f}")
            with col_rt3:
                # Count restriction sites
                enzymes = enzyme_pair.split(" / ")
                total_sites = 0
                for enzyme in enzymes:
                    sites = gsynth_find_restriction_sites(clean_vector, enzyme)
                    total_sites += len(sites)
                st.metric("RE Sites", total_sites)
    
    with col_input2:
        st.markdown("#### Insert Sequences")
        
        # ✅ FIXED: Changed height from 60 to 70 (minimum requirement)
        forward_sequence = st.text_area(
            "Forward strand (5'→3')",
            value=st.session_state.get("example_forward", ""),
            height=70,  # ✅ Fixed: was 60, now 70 (meets minimum 68px requirement)
            placeholder="Forward DNA strand...",
            help="Enter the forward insert strand with sticky ends"
        )
        
        # ✅ FIXED: Changed height from 60 to 70 (minimum requirement)
        reverse_sequence = st.text_area(
            "Reverse strand (5'→3')",
            value=st.session_state.get("example_reverse", ""),
            height=70,  # ✅ Fixed: was 60, now 70 (meets minimum 68px requirement)
            placeholder="Reverse DNA strand...",
            help="Enter the reverse insert strand with sticky ends"
        )
    
    # Action buttons
    st.markdown("#### 🧪 Analysis Controls")
    col_btn1, col_btn2 = st.columns(2)
    
    with col_btn1:
        if st.button("🔬 Digest Vector", type="primary", use_container_width=True):
            if not vector_sequence.strip():
                st.error("❌ Please enter a vector sequence")
            else:
                clean_vector = clean_dna_sequence(vector_sequence) if BIO_UTILS_AVAILABLE else vector_sequence.upper()
                
                with st.spinner("🧬 Performing restriction digestion..."):
                    result = gsynth_simulate_digestion_enhanced(clean_vector, enzyme_pair)
                    
                    if "error" in result:
                        st.error(f"❌ Digestion failed: {result['error']}")
                    else:
                        st.session_state["ligation_data"]["digestion_result"] = result
                        st.session_state["ligation_data"]["plasmid_name"] = plasmid_name
                        st.session_state["ligation_data"]["insert_name"] = insert_name
                        st.success("✅ Vector digestion completed!")
                        
                        fragments = result["fragments"]
                        st.info(f"""
                        **Digestion Results:**
                        - Original: {len(result['original_vector']):,} bp
                        - Linearized: {fragments['vector_length']:,} bp
                        - Fragment removed: {fragments['fragment_length']:,} bp
                        - Recovery: {(fragments['total_recovered']/len(result['original_vector'])*100):.1f}%
                        """)
    
    with col_btn2:
        if st.button("🔗 Check Ligation", type="secondary", use_container_width=True):
            if not st.session_state["ligation_data"].get("digestion_result"):
                st.error("❌ Please perform vector digestion first")
            elif not forward_sequence.strip() or not reverse_sequence.strip():
                st.error("❌ Please enter both forward and reverse sequences")
            else:
                clean_forward = clean_dna_sequence(forward_sequence) if BIO_UTILS_AVAILABLE else forward_sequence.upper()
                clean_reverse = clean_dna_sequence(reverse_sequence) if BIO_UTILS_AVAILABLE else reverse_sequence.upper()
                
                with st.spinner("🧬 Analyzing ligation compatibility..."):
                    sticky_result = gsynth_validate_sticky_ends_enhanced(clean_forward, clean_reverse, enzyme_pair)
                    hybrid_result = gsynth_check_hybridization_enhanced(clean_forward, clean_reverse, enzyme_pair)
                    
                    validation_result = {
                        **sticky_result,
                        **hybrid_result,
                        "enzyme_pair": enzyme_pair,
                        "forward_sequence": clean_forward,
                        "reverse_sequence": clean_reverse
                    }
                    
                    st.session_state["ligation_data"]["validation_result"] = validation_result
                    
                    if validation_result.get("overall_valid") and validation_result.get("hybridization_success"):
                        st.success("🎉 Ligation compatibility confirmed!")
                        st.balloons()
                    else:
                        st.error("❌ Ligation compatibility issues detected")
    
    st.markdown("</div>", unsafe_allow_html=True)
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # SECTION 2: VISUALIZATION RESULTS SECTION
    # ═══════════════════════════════════════════════════════════════════════════════
    
    # Show results only if we have data
    digestion_result = st.session_state["ligation_data"].get("digestion_result")
    validation_result = st.session_state["ligation_data"].get("validation_result")
    stored_plasmid_name = st.session_state["ligation_data"].get("plasmid_name", plasmid_name)
    stored_insert_name = st.session_state["ligation_data"].get("insert_name", insert_name)
    
    if digestion_result:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, #F0FDF4, #DCFCE7);
            padding: 2rem;
            border-radius: 20px;
            margin-bottom: 2rem;
            border: 3px solid #22C55E;
            box-shadow: 0 8px 25px rgba(34, 197, 94, 0.2);
        ">
            <h2 style="color: #15803D; margin: 0 0 1.5rem 0; font-size: 1.8rem;">
                🔬 Visualization Results
            </h2>
        """, unsafe_allow_html=True)
        
        # Status indicators
        col_status1, col_status2, col_status3 = st.columns(3)
        
        with col_status1:
            st.metric("Vector Status", "✅ Digested")
        
        with col_status2:
            if validation_result:
                overall_success = validation_result.get("overall_valid", False) and validation_result.get("hybridization_success", False)
                st.metric("Ligation Status", "✅ Compatible" if overall_success else "❌ Issues")
            else:
                st.metric("Ligation Status", "⏳ Pending")
        
        with col_status3:
            if validation_result and validation_result.get("overall_valid"):
                insert_size = len(validation_result.get("forward_sequence", ""))
                st.metric("Insert Size", f"{insert_size:,} bp")
            else:
                st.metric("Insert Size", "N/A")
        
        # PLASMID VISUALIZATION (Moved here as requested)
        st.markdown("#### 🧬 Plasmid Map Visualization")
        
        # Determine success for visualization
        ligation_success = False
        if validation_result:
            ligation_success = (validation_result.get("overall_valid", False) and 
                              validation_result.get("hybridization_success", False))
        
        # Create enhanced plasmid map
        plasmid_fig = create_enhanced_plasmid_map(
            digestion_result["linearized_vector"],
            digestion_result["cut_positions"],
            enzyme_pair,
            validation_result["forward_sequence"] if ligation_success and validation_result else None,
            ligation_success,
            stored_plasmid_name,
            stored_insert_name
        )
        
        st.plotly_chart(plasmid_fig, use_container_width=True)
        
        # Additional visualizations in tabs
        if validation_result:
            st.markdown("#### 📊 Detailed Analysis")
            
            tab1, tab2 = st.tabs([
                "🔗 Junction Analysis", 
                "📋 Analysis Report"
            ])
            
            with tab1:
                junction_fig = create_enhanced_junction_analysis(
                    validation_result["forward_sequence"],
                    validation_result["reverse_sequence"],
                    enzyme_pair,
                    validation_result,
                    stored_insert_name
                )
                
                st.plotly_chart(junction_fig, use_container_width=True)
                
                # Detailed analysis results
                if validation_result.get("forward_analysis"):
                    fwd_analysis = validation_result["forward_analysis"]
                    rev_analysis = validation_result["reverse_analysis"]
                    
                    col_detail1, col_detail2 = st.columns(2)
                    
                    with col_detail1:
                        st.markdown(f"""
                        **Forward End Analysis ({stored_insert_name}):**
                        - Expected: `{fwd_analysis['expected']}`
                        - Found: `{fwd_analysis['found']}`
                        - Similarity: {fwd_analysis['similarity']:.1f}%
                        - Status: {'✅ Valid' if fwd_analysis['valid'] else '❌ Invalid'}
                        """)
                    
                    with col_detail2:
                        st.markdown(f"""
                        **Reverse End Analysis ({stored_insert_name}):**
                        - Expected: `{rev_analysis['expected']}`
                        - Found: `{rev_analysis['found']}`
                        - Similarity: {rev_analysis['similarity']:.1f}%
                        - Status: {'✅ Valid' if rev_analysis['valid'] else '❌ Invalid'}
                        """)
            
            with tab2:
                # Generate comprehensive report
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                ligation_success = (validation_result.get("overall_valid", False) and 
                                  validation_result.get("hybridization_success", False))
                
                report = f"""
# SYNTHETIC INSERT LIGATION ANALYSIS REPORT
==========================================

**Analysis Date:** {timestamp}
**Plasmid Name:** {stored_plasmid_name}
**Insert Name:** {stored_insert_name}
**Enzyme Pair:** {enzyme_pair}
**Analysis Status:** {'✅ SUCCESSFUL LIGATION' if ligation_success else '❌ LIGATION ISSUES'}

## MOLECULAR SPECIFICATIONS
=========================
**Plasmid Vector:**
- Name: {stored_plasmid_name}
- Original Length: {len(digestion_result['original_vector']):,} bp
- Linearized Length: {digestion_result['fragments']['vector_length']:,} bp
- GC Content: {calculate_gc(digestion_result['linearized_vector']) if BIO_UTILS_AVAILABLE else 'N/A':.1f}%

**Insert Fragment:**
- Name: {stored_insert_name}
- Length: {len(validation_result['forward_sequence']):,} bp
- GC Content: {calculate_gc(validation_result['forward_sequence']) if BIO_UTILS_AVAILABLE else 'N/A':.1f}%

**Final Recombinant:**
- Expected Size: {digestion_result['fragments']['vector_length'] + len(validation_result['forward_sequence']):,} bp
- Insert Integration: {'Successful' if ligation_success else 'Failed'}

## RESTRICTION ENZYME ANALYSIS
============================
**Enzyme Pair:** {enzyme_pair}
**Recognition Sites Found:** {len(digestion_result['forward_sites']) + len(digestion_result['reverse_sites'])}
**Cut Positions:** {digestion_result['cut_positions']}

## LIGATION COMPATIBILITY
=======================
**Sticky End Validation:** {'✅ PASSED' if validation_result.get('overall_valid') else '❌ FAILED'}
**Hybridization Analysis:** {'✅ PASSED' if validation_result.get('hybridization_success') else '❌ FAILED'}
**Match Percentage:** {validation_result.get('match_percentage', 0):.1f}%

## EXPERIMENTAL RECOMMENDATION
============================
**Protocol:** {'Proceed with standard ligation protocol' if ligation_success else 'Optimize sequences before attempting ligation'}
**Conditions:** T4 DNA Ligase, 16°C overnight or room temperature for 2-4 hours
**Expected Success Rate:** {'High (>85%)' if ligation_success else 'Low (<25%)'}

Generated by Synthetic Insert Ligation Check v2.0
Enhanced Layout Version with Manual Name Input and Fixed Text Area Heights
                """
                
                st.text_area("📄 Complete Analysis Report", report, height=400)
                
                # Download button
                st.download_button(
                    "📄 Download Complete Report",
                    report,
                    f"ligation_report_{stored_plasmid_name}_{stored_insert_name}_{timestamp.replace(' ', '_').replace(':', '-')}.txt",
                    "text/plain",
                    use_container_width=True
                )
        
        st.markdown("</div>", unsafe_allow_html=True)
    
    else:
        # Show placeholder when no analysis has been run
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, #FFFBEB, #FEF3C7);
            padding: 2rem;
            border-radius: 20px;
            margin-bottom: 2rem;
            border: 3px solid #F59E0B;
            box-shadow: 0 8px 25px rgba(245, 158, 11, 0.2);
            text-align: center;
        ">
            <h2 style="color: #92400E; margin: 0 0 1rem 0; font-size: 1.8rem;">
                🔬 Visualization Results
            </h2>
            <p style="color: #92400E; font-size: 1.1rem; margin: 0;">
                Run vector digestion and ligation check to see visualization results here
            </p>
        </div>
        """, unsafe_allow_html=True)

def app():
    """Entry point for modular integration"""
    main()

if __name__ == "__main__":
    main()