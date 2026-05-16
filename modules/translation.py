# modules/translation.py

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from collections import Counter
import io
import re

# ───────────────────────────────────────────────────────────────────────────────
# GENETIC CODE AND CODON TABLES
# ───────────────────────────────────────────────────────────────────────────────

# Standard genetic code
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Amino acid properties
AA_PROPERTIES = {
    'A': {'name': 'Alanine', 'mw': 89.1, 'hydrophobic': True, 'charged': False},
    'R': {'name': 'Arginine', 'mw': 174.2, 'hydrophobic': False, 'charged': True},
    'N': {'name': 'Asparagine', 'mw': 132.1, 'hydrophobic': False, 'charged': False},
    'D': {'name': 'Aspartic acid', 'mw': 133.1, 'hydrophobic': False, 'charged': True},
    'C': {'name': 'Cysteine', 'mw': 121.2, 'hydrophobic': False, 'charged': False},
    'E': {'name': 'Glutamic acid', 'mw': 147.1, 'hydrophobic': False, 'charged': True},
    'Q': {'name': 'Glutamine', 'mw': 146.1, 'hydrophobic': False, 'charged': False},
    'G': {'name': 'Glycine', 'mw': 75.1, 'hydrophobic': False, 'charged': False},
    'H': {'name': 'Histidine', 'mw': 155.2, 'hydrophobic': False, 'charged': True},
    'I': {'name': 'Isoleucine', 'mw': 131.2, 'hydrophobic': True, 'charged': False},
    'L': {'name': 'Leucine', 'mw': 131.2, 'hydrophobic': True, 'charged': False},
    'K': {'name': 'Lysine', 'mw': 146.2, 'hydrophobic': False, 'charged': True},
    'M': {'name': 'Methionine', 'mw': 149.2, 'hydrophobic': True, 'charged': False},
    'F': {'name': 'Phenylalanine', 'mw': 165.2, 'hydrophobic': True, 'charged': False},
    'P': {'name': 'Proline', 'mw': 115.1, 'hydrophobic': False, 'charged': False},
    'S': {'name': 'Serine', 'mw': 105.1, 'hydrophobic': False, 'charged': False},
    'T': {'name': 'Threonine', 'mw': 119.1, 'hydrophobic': False, 'charged': False},
    'W': {'name': 'Tryptophan', 'mw': 204.2, 'hydrophobic': True, 'charged': False},
    'Y': {'name': 'Tyrosine', 'mw': 181.2, 'hydrophobic': False, 'charged': False},
    'V': {'name': 'Valine', 'mw': 117.1, 'hydrophobic': True, 'charged': False},
    '*': {'name': 'Stop', 'mw': 0, 'hydrophobic': False, 'charged': False}
}

# Three-letter to one-letter conversion
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'STOP': '*', 'TER': '*', 'END': '*'
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items() if k != 'STOP' and k != 'TER' and k != 'END'}
ONE_TO_THREE['*'] = 'STOP'

# Organism-specific codon usage tables
CODON_USAGE_TABLES = {
    'E. coli': {
        'F': {'TTT': 0.58, 'TTC': 0.42},
        'L': {'TTG': 0.13, 'TTA': 0.14, 'CTG': 0.50, 'CTA': 0.04, 'CTC': 0.10, 'CTT': 0.10},
        'S': {'TCG': 0.15, 'TCA': 0.14, 'TCC': 0.15, 'TCT': 0.15, 'AGC': 0.28, 'AGT': 0.15},
        'Y': {'TAT': 0.59, 'TAC': 0.41},
        'C': {'TGT': 0.46, 'TGC': 0.54},
        'W': {'TGG': 1.00},
        'P': {'CCG': 0.53, 'CCA': 0.20, 'CCC': 0.12, 'CCT': 0.16},
        'H': {'CAT': 0.57, 'CAC': 0.43},
        'Q': {'CAG': 0.65, 'CAA': 0.35},
        'R': {'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.11, 'AGA': 0.07, 'AGG': 0.04},
        'I': {'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11},
        'M': {'ATG': 1.00},
        'T': {'ACG': 0.27, 'ACA': 0.17, 'ACC': 0.40, 'ACT': 0.16},
        'N': {'AAT': 0.49, 'AAC': 0.51},
        'K': {'AAA': 0.74, 'AAG': 0.26},
        'V': {'GTG': 0.37, 'GTA': 0.15, 'GTC': 0.20, 'GTT': 0.28},
        'A': {'GCG': 0.36, 'GCA': 0.21, 'GCC': 0.26, 'GCT': 0.18},
        'D': {'GAT': 0.63, 'GAC': 0.37},
        'E': {'GAA': 0.68, 'GAG': 0.32},
        'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15},
        '*': {'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30}
    },
    'Human': {
        'F': {'TTT': 0.46, 'TTC': 0.54},
        'L': {'TTG': 0.13, 'TTA': 0.08, 'CTG': 0.40, 'CTA': 0.07, 'CTC': 0.20, 'CTT': 0.13},
        'S': {'TCG': 0.05, 'TCA': 0.12, 'TCC': 0.22, 'TCT': 0.19, 'AGC': 0.24, 'AGT': 0.15},
        'Y': {'TAT': 0.44, 'TAC': 0.56},
        'C': {'TGT': 0.46, 'TGC': 0.54},
        'W': {'TGG': 1.00},
        'P': {'CCG': 0.07, 'CCA': 0.27, 'CCC': 0.32, 'CCT': 0.29},
        'H': {'CAT': 0.42, 'CAC': 0.58},
        'Q': {'CAG': 0.75, 'CAA': 0.25},
        'R': {'CGT': 0.08, 'CGC': 0.18, 'CGA': 0.11, 'CGG': 0.20, 'AGA': 0.20, 'AGG': 0.20},
        'I': {'ATT': 0.36, 'ATC': 0.47, 'ATA': 0.17},
        'M': {'ATG': 1.00},
        'T': {'ACG': 0.12, 'ACA': 0.24, 'ACC': 0.36, 'ACT': 0.25},
        'N': {'AAT': 0.47, 'AAC': 0.53},
        'K': {'AAA': 0.43, 'AAG': 0.57},
        'V': {'GTG': 0.47, 'GTA': 0.11, 'GTC': 0.24, 'GTT': 0.18},
        'A': {'GCG': 0.11, 'GCA': 0.23, 'GCC': 0.40, 'GCT': 0.27},
        'D': {'GAT': 0.46, 'GAC': 0.54},
        'E': {'GAA': 0.42, 'GAG': 0.58},
        'G': {'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25},
        '*': {'TAA': 0.30, 'TAG': 0.24, 'TGA': 0.46}
    },
    'S. cerevisiae': {
        'F': {'TTT': 0.59, 'TTC': 0.41},
        'L': {'TTG': 0.29, 'TTA': 0.28, 'CTG': 0.11, 'CTA': 0.14, 'CTC': 0.06, 'CTT': 0.13},
        'S': {'TCG': 0.10, 'TCA': 0.21, 'TCC': 0.16, 'TCT': 0.26, 'AGC': 0.11, 'AGT': 0.16},
        'Y': {'TAT': 0.56, 'TAC': 0.44},
        'C': {'TGT': 0.63, 'TGC': 0.37},
        'W': {'TGG': 1.00},
        'P': {'CCG': 0.12, 'CCA': 0.42, 'CCC': 0.15, 'CCT': 0.31},
        'H': {'CAT': 0.64, 'CAC': 0.36},
        'Q': {'CAG': 0.31, 'CAA': 0.69},
        'R': {'CGT': 0.15, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.48, 'AGG': 0.21},
        'I': {'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27},
        'M': {'ATG': 1.00},
        'T': {'ACG': 0.14, 'ACA': 0.35, 'ACC': 0.22, 'ACT': 0.29},
        'N': {'AAT': 0.59, 'AAC': 0.41},
        'K': {'AAA': 0.58, 'AAG': 0.42},
        'V': {'GTG': 0.12, 'GTA': 0.21, 'GTC': 0.15, 'GTT': 0.52},
        'A': {'GCG': 0.11, 'GCA': 0.29, 'GCC': 0.22, 'GCT': 0.38},
        'D': {'GAT': 0.65, 'GAC': 0.35},
        'E': {'GAA': 0.70, 'GAG': 0.30},
        'G': {'GGT': 0.47, 'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12},
        '*': {'TAA': 0.47, 'TAG': 0.23, 'TGA': 0.30}
    }
}

# ───────────────────────────────────────────────────────────────────────────────
# ENHANCED CSS FOR TRANSLATION MODULE
# ───────────────────────────────────────────────────────────────────────────────

def inject_translation_css():
    """Enhanced CSS specifically for Translation module with red accent theme"""
    st.markdown("""
    <style>
    /* Translation-specific enhancements */
    .translation-input-card {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%) !important;
        border: 2px solid #e2e8f0 !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08) !important;
        transition: all 0.3s ease !important;
    }
    
    .translation-input-card:hover {
        border-color: #ff4b4b !important;
        box-shadow: 0 8px 25px rgba(255, 75, 75, 0.15) !important;
        transform: translateY(-2px) !important;
    }
    
    .translation-sequence-display {
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
    
    .translation-results-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 1rem 0 !important;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.3) !important;
        text-align: center !important;
    }
    
    .translation-section-header {
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
    
    .translation-metric-container {
        background: rgba(255, 255, 255, 0.2) !important;
        border-radius: 12px !important;
        padding: 1rem !important;
        margin: 0.5rem 0 !important;
        border: 1px solid rgba(255, 255, 255, 0.3) !important;
        backdrop-filter: blur(10px) !important;
    }
    
    .translation-info-box {
        background: linear-gradient(135deg, rgba(34, 197, 94, 0.1), rgba(16, 185, 129, 0.1)) !important;
        border-left: 4px solid #10b981 !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #059669 !important;
    }
    
    .translation-warning-box {
        background: linear-gradient(135deg, rgba(251, 146, 60, 0.1), rgba(245, 158, 11, 0.1)) !important;
        border-left: 4px solid #f59e0b !important;
        border-radius: 8px !important;
        padding: 1rem !important;
        margin: 1rem 0 !important;
        color: #d97706 !important;
    }
    
    .translation-copy-button {
        background: linear-gradient(135deg, #ff4b4b, #cc0000) !important;
        color: white !important;
        border: none !important;
        padding: 8px 16px !important;
        border-radius: 8px !important;
        cursor: pointer !important;
        font-weight: 600 !important;
        font-size: 14px !important;
        transition: all 0.3s ease !important;
        margin: 0.5rem 0 !important;
    }
    
    .translation-copy-button:hover {
        transform: scale(1.05) !important;
        box-shadow: 0 4px 15px rgba(255, 75, 75, 0.3) !important;
    }
    </style>
    """, unsafe_allow_html=True)

# ───────────────────────────────────────────────────────────────────────────────
# TRANSLATION FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def clean_dna_sequence(seq):
    """Clean and validate DNA sequence"""
    if not seq:
        return ""
    
    # Remove whitespace and convert to uppercase
    clean_seq = re.sub(r'\s+', '', seq.upper())
    
    # Remove non-DNA characters
    clean_seq = re.sub(r'[^ATGCNRYSWKMBDHV]', '', clean_seq)
    
    return clean_seq

def translate_sequence(dna_seq, frame=0, find_start=True):
    """Enhanced translation with multiple options"""
    if not dna_seq:
        return ""
    
    # Clean sequence
    clean_seq = clean_dna_sequence(dna_seq)
    
    # Adjust for reading frame
    seq_to_translate = clean_seq[frame:]
    
    # Find start codon if requested
    if find_start:
        start_pos = seq_to_translate.find('ATG')
        if start_pos != -1:
            seq_to_translate = seq_to_translate[start_pos:]
    
    # Translate
    protein = ""
    for i in range(0, len(seq_to_translate), 3):
        codon = seq_to_translate[i:i+3]
        if len(codon) == 3:
            aa = GENETIC_CODE.get(codon, 'X')
            protein += aa
            if aa == '*':  # Stop at first stop codon
                break
    
    return protein

def translate_all_frames(dna_seq):
    """Translate all 6 reading frames"""
    frames = {}
    
    # Forward frames
    for frame in range(3):
        frames[f'Frame +{frame+1}'] = translate_sequence(dna_seq, frame, find_start=False)
    
    # Reverse complement frames
    reverse_seq = str(reverse_complement(dna_seq))
    for frame in range(3):
        frames[f'Frame -{frame+1}'] = translate_sequence(reverse_seq, frame, find_start=False)
    
    return frames

def reverse_complement(dna_seq):
    """Calculate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(dna_seq))

def reverse_translate_to_dna(protein_seq, organism='E. coli', method='most_frequent'):
    """Enhanced reverse translation with organism-specific codon usage"""
    if not protein_seq:
        return ""
    
    codon_table = CODON_USAGE_TABLES.get(organism, CODON_USAGE_TABLES['E. coli'])
    dna_seq = ""
    
    for aa in protein_seq.upper():
        if aa in codon_table:
            codons = codon_table[aa]
            
            if method == 'most_frequent':
                # Choose most frequent codon
                best_codon = max(codons.items(), key=lambda x: x[1])[0]
            elif method == 'random':
                # Choose random codon
                best_codon = np.random.choice(list(codons.keys()))
            elif method == 'weighted':
                # Choose based on frequency distribution
                codon_list = list(codons.keys())
                frequencies = list(codons.values())
                best_codon = np.random.choice(codon_list, p=frequencies)
            else:
                best_codon = list(codons.keys())[0]
            
            dna_seq += best_codon
        elif aa == '*':
            # Handle stop codons
            stop_codons = codon_table.get('*', {'TAA': 1.0})
            best_codon = max(stop_codons.items(), key=lambda x: x[1])[0]
            dna_seq += best_codon
    
    return dna_seq

def convert_to_three_letter(one_letter_seq):
    """Convert one-letter amino acid code to three-letter"""
    if not one_letter_seq:
        return ""
    
    three_letter = []
    for aa in one_letter_seq.upper():
        if aa in ONE_TO_THREE:
            three_letter.append(ONE_TO_THREE[aa])
        else:
            three_letter.append('XXX')
    
    return ' '.join(three_letter)

def convert_three_to_one(three_letter_seq):
    """Convert three-letter amino acid code to one-letter"""
    if not three_letter_seq:
        return ""
    
    # Split by spaces and clean
    three_codes = [code.strip().upper() for code in three_letter_seq.split()]
    
    one_letter = ""
    for code in three_codes:
        if code in THREE_TO_ONE:
            one_letter += THREE_TO_ONE[code]
        else:
            one_letter += 'X'
    
    return one_letter

def calculate_gc_content(dna_seq):
    """Calculate GC content percentage"""
    if not dna_seq:
        return 0
    
    gc_count = dna_seq.count('G') + dna_seq.count('C')
    return (gc_count / len(dna_seq)) * 100

def analyze_protein_properties(protein_seq):
    """Analyze protein sequence properties"""
    if not protein_seq:
        return {}
    
    properties = {
        'length': len(protein_seq),
        'molecular_weight': 0,
        'hydrophobic_residues': 0,
        'charged_residues': 0,
        'composition': Counter(protein_seq)
    }
    
    for aa in protein_seq:
        if aa in AA_PROPERTIES:
            properties['molecular_weight'] += AA_PROPERTIES[aa]['mw']
            if AA_PROPERTIES[aa]['hydrophobic']:
                properties['hydrophobic_residues'] += 1
            if AA_PROPERTIES[aa]['charged']:
                properties['charged_residues'] += 1
    
    # Calculate percentages
    if properties['length'] > 0:
        properties['hydrophobic_percent'] = (properties['hydrophobic_residues'] / properties['length']) * 100
        properties['charged_percent'] = (properties['charged_residues'] / properties['length']) * 100
    
    return properties

# ───────────────────────────────────────────────────────────────────────────────
# ENHANCED VISUALIZATION FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def create_reading_frame_visualization(frames_data):
    """Create interactive reading frame visualization with red theme"""
    
    fig = make_subplots(
        rows=6, cols=1,
        subplot_titles=list(frames_data.keys()),
        vertical_spacing=0.02,
        shared_xaxes=True
    )
    
    colors = ['#ff4b4b', '#cc0000', '#ff6b6b', '#e53e3e', '#ff8a80', '#ffcdd2']
    
    for i, (frame_name, protein_seq) in enumerate(frames_data.items()):
        # Create sequence data for visualization
        x_positions = list(range(len(protein_seq)))
        y_values = [1] * len(protein_seq)
        
        # Color code amino acids
        aa_colors = []
        for aa in protein_seq:
            if aa == '*':
                aa_colors.append('#FF1744')  # Red for stop codons
            elif aa == 'M':
                aa_colors.append('#4CAF50')  # Green for methionine
            elif aa in 'ACFGILMPVW':  # Hydrophobic
                aa_colors.append('#ff4b4b')  # Red theme
            elif aa in 'DEHKR':  # Charged
                aa_colors.append('#FF9800')  # Orange
            else:
                aa_colors.append('#9E9E9E')  # Gray for others
        
        fig.add_trace(
            go.Scatter(
                x=x_positions,
                y=y_values,
                mode='markers',
                marker=dict(
                    color=aa_colors,
                    size=8,
                    line=dict(width=1, color='white')
                ),
                text=[f"Position {i+1}: {aa}" for i, aa in enumerate(protein_seq)],
                hovertemplate='%{text}<extra></extra>',
                name=frame_name,
                showlegend=False
            ),
            row=i+1, col=1
        )
    
    fig.update_layout(
        height=800,
        title="Interactive Reading Frame Analysis",
        plot_bgcolor='#f8f9ff',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12)
    )
    
    # Update y-axes
    for i in range(6):
        fig.update_yaxes(showticklabels=False, row=i+1, col=1)
    
    fig.update_xaxes(title_text="Amino Acid Position", row=6, col=1)
    
    return fig

def create_amino_acid_composition_chart(protein_seq):
    """Create amino acid composition chart with red theme"""
    
    composition = Counter(protein_seq)
    
    # Prepare data
    amino_acids = list(composition.keys())
    counts = list(composition.values())
    percentages = [count/len(protein_seq)*100 for count in counts]
    
    # Create color map based on properties with red theme
    colors = []
    for aa in amino_acids:
        if aa == '*':
            colors.append('#FF1744')
        elif aa in 'ACFGILMPVW':  # Hydrophobic
            colors.append('#ff4b4b')
        elif aa in 'DEHKR':  # Charged
            colors.append('#FF9800')
        elif aa in 'NQST':  # Polar
            colors.append('#4CAF50')
        else:
            colors.append('#9E9E9E')
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=amino_acids,
        y=percentages,
        marker_color=colors,
        text=[f"{p:.1f}%" for p in percentages],
        textposition='auto',
        hovertemplate='<b>%{x}</b><br>Count: %{customdata}<br>Percentage: %{y:.1f}%<extra></extra>',
        customdata=counts
    ))
    
    fig.update_layout(
        title="Amino Acid Composition Analysis",
        xaxis_title="Amino Acid",
        yaxis_title="Percentage (%)",
        plot_bgcolor='#f8f9ff',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12),
        height=400
    )
    
    return fig

def create_copy_button_html(sequence, label, button_id):
    """Create enhanced copy button with red theme"""
    return f"""
    <div style="margin: 10px 0;">
        <button onclick="copyToClipboard{button_id}()" class="translation-copy-button">
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

# ───────────────────────────────────────────────────────────────────────────────
# MAIN APPLICATION
# ───────────────────────────────────────────────────────────────────────────────

def main():
    """Enhanced Translation & Reverse Translation tool with red theme"""
    
    # Inject enhanced CSS
    inject_translation_css()
    
    # Professional header with red theme (matching your design preference)
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
            🧬 Translation & Reverse Translation
        </h2>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.1rem; opacity: 0.9;">
            Professional DNA-protein translation tools with codon optimization and analysis
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Enhanced tabs with modern design
    tab1, tab2, tab3, tab4 = st.tabs([
        "🧬 DNA → Protein", 
        "🔄 Protein → DNA", 
        "🔤 AA Conversion", 
        "📊 Analysis"
    ])
    
    with tab1:
        st.markdown('<div class="translation-section-header">🧬 DNA to Protein Translation</div>', unsafe_allow_html=True)
        
        st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            dna_input = st.text_area(
                "**Enter DNA Sequence (5' → 3')**",
                height=150,
                placeholder="ATGAAGTGGGTAACCTTCATC...",
                help="Enter DNA nucleotide sequence (A, T, G, C)"
            )
        
        with col2:
            reading_frame_options = ["None (Start from beginning)"] + [f"Frame +{i+1}" for i in range(3)]
            reading_frame_selection = st.selectbox(
                "**Reading Frame**",
                reading_frame_options,
                index=0,
                help="Choose the reading frame for translation or None to start from beginning"
            )
            
            # Extract frame number or set to None
            if reading_frame_selection == "None (Start from beginning)":
                reading_frame = 0
                find_start_codon = False
            else:
                reading_frame = int(reading_frame_selection.split("+")[1]) - 1
                find_start_codon = st.checkbox(
                    "**Find first ATG start codon**",
                    value=True,
                    help="Start translation from the first ATG"
                )
            
            show_all_frames = st.checkbox(
                "**Show all 6 reading frames**",
                value=False,
                help="Display translation in all possible reading frames"
            )
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Process button
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            translate_clicked = st.button("🚀 Translate DNA", type="primary", use_container_width=True)
        
        if translate_clicked:
            if not dna_input.strip():
                st.error("Please enter a DNA sequence.")
            else:
                clean_dna = clean_dna_sequence(dna_input)
                if not clean_dna:
                    st.error("No valid DNA characters found.")
                else:
                    with st.spinner("Translating sequence..."):
                        st.success("✅ Translation completed successfully!")
                        
                        if show_all_frames:
                            # Show all 6 reading frames
                            st.markdown('<div class="translation-section-header">📊 All Reading Frames Analysis</div>', unsafe_allow_html=True)
                            
                            all_frames = translate_all_frames(clean_dna)
                            
                            # Create interactive visualization
                            frames_fig = create_reading_frame_visualization(all_frames)
                            st.plotly_chart(frames_fig, use_container_width=True)
                            
                            # Display results in table
                            frames_df = pd.DataFrame([
                                {
                                    'Frame': frame,
                                    'Protein Sequence': protein[:50] + '...' if len(protein) > 50 else protein,
                                    'Length': len(protein),
                                    'Stops': protein.count('*')
                                }
                                for frame, protein in all_frames.items()
                            ])
                            
                            st.dataframe(frames_df, use_container_width=True)
                            
                            # Copy buttons for each frame
                            st.markdown('<div class="translation-section-header">📋 Copy Individual Frames</div>', unsafe_allow_html=True)
                            cols = st.columns(3)
                            for i, (frame, protein) in enumerate(all_frames.items()):
                                with cols[i % 3]:
                                    st.markdown(f"**{frame}**")
                                    copy_html = create_copy_button_html(protein, frame, f"frame{i}")
                                    st.components.v1.html(copy_html, height=80)
                        
                        else:
                            # Single frame translation
                            protein_seq = translate_sequence(clean_dna, reading_frame, find_start_codon)
                            
                            st.markdown('<div class="translation-section-header">🧬 Translation Result</div>', unsafe_allow_html=True)
                            
                            col_result1, col_result2 = st.columns(2)
                            
                            with col_result1:
                                st.markdown("**DNA Sequence**")
                                st.markdown(f'<div class="translation-sequence-display">{clean_dna}</div>', unsafe_allow_html=True)
                                
                                copy_dna_html = create_copy_button_html(clean_dna, "DNA", "dna1")
                                st.components.v1.html(copy_dna_html, height=80)
                            
                            with col_result2:
                                st.markdown("**Protein Sequence**")
                                st.markdown(f'<div class="translation-sequence-display">{protein_seq}</div>', unsafe_allow_html=True)
                                
                                copy_protein_html = create_copy_button_html(protein_seq, "Protein", "protein1")
                                st.components.v1.html(copy_protein_html, height=80)
                            
                            # Analysis metrics
                            st.markdown('<div class="translation-section-header">📊 Translation Metrics</div>', unsafe_allow_html=True)
                            
                            st.markdown('<div class="translation-results-card">', unsafe_allow_html=True)
                            
                            col_metric1, col_metric2, col_metric3, col_metric4 = st.columns(4)
                            
                            with col_metric1:
                                st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                st.metric("DNA Length", f"{len(clean_dna)} bp")
                                st.markdown('</div>', unsafe_allow_html=True)
                            
                            with col_metric2:
                                st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                st.metric("Protein Length", f"{len(protein_seq)} aa")
                                st.markdown('</div>', unsafe_allow_html=True)
                            
                            with col_metric3:
                                st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                st.metric("GC Content", f"{calculate_gc_content(clean_dna):.1f}%")
                                st.markdown('</div>', unsafe_allow_html=True)
                            
                            with col_metric4:
                                st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                stop_count = protein_seq.count('*')
                                st.metric("Stop Codons", stop_count)
                                st.markdown('</div>', unsafe_allow_html=True)
                            
                            st.markdown('</div>', unsafe_allow_html=True)
                            
                            # Protein properties analysis
                            if protein_seq:
                                properties = analyze_protein_properties(protein_seq)
                                
                                st.markdown('<div class="translation-section-header">🔬 Protein Properties</div>', unsafe_allow_html=True)
                                
                                st.markdown('<div class="translation-results-card">', unsafe_allow_html=True)
                                prop_col1, prop_col2, prop_col3 = st.columns(3)
                                
                                with prop_col1:
                                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                    st.metric("Molecular Weight", f"{properties['molecular_weight']:.1f} Da")
                                    st.markdown('</div>', unsafe_allow_html=True)
                                
                                with prop_col2:
                                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                    st.metric("Hydrophobic %", f"{properties.get('hydrophobic_percent', 0):.1f}%")
                                    st.markdown('</div>', unsafe_allow_html=True)
                                
                                with prop_col3:
                                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                                    st.metric("Charged %", f"{properties.get('charged_percent', 0):.1f}%")
                                    st.markdown('</div>', unsafe_allow_html=True)
                                
                                st.markdown('</div>', unsafe_allow_html=True)
                                
                                # Amino acid composition chart
                                if len(protein_seq) > 1:
                                    composition_fig = create_amino_acid_composition_chart(protein_seq)
                                    st.plotly_chart(composition_fig, use_container_width=True)
    
    with tab2:
        st.markdown('<div class="translation-section-header">🔄 Protein to DNA Reverse Translation</div>', unsafe_allow_html=True)
        
        st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            protein_input = st.text_area(
                "**Enter Protein Sequence (One-letter codes)**",
                height=150,
                placeholder="MKWVTFISLLLLFSSAYSRGVFRRD...",
                help="Enter protein sequence using single-letter amino acid codes"
            )
        
        with col2:
            organism_options = ["None (No optimization)"] + list(CODON_USAGE_TABLES.keys())
            organism_selection = st.selectbox(
                "**Target Organism**",
                organism_options,
                index=1,  # Default to E. coli
                help="Choose organism for codon optimization or None for no optimization"
            )
            
            organism = None if organism_selection == "None (No optimization)" else organism_selection
            
            if organism:
                reverse_method = st.selectbox(
                    "**Codon Selection Method**",
                    ['most_frequent', 'weighted', 'random'],
                    format_func=lambda x: {
                        'most_frequent': 'Most Frequent Codon',
                        'weighted': 'Frequency-Weighted Random',
                        'random': 'Random Selection'
                    }[x],
                    help="Method for choosing codons during reverse translation"
                )
            else:
                reverse_method = 'most_frequent'
                st.markdown('<div class="translation-info-box">ℹ️ No organism selected - using standard genetic code without optimization</div>', unsafe_allow_html=True)
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Process button
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            reverse_translate_clicked = st.button("🚀 Reverse Translate", type="primary", use_container_width=True)
        
        if reverse_translate_clicked:
            if not protein_input.strip():
                st.error("Please enter a protein sequence.")
            else:
                clean_protein = protein_input.strip().upper()
                
                # Validate protein sequence
                invalid_chars = set(clean_protein) - set('ACDEFGHIKLMNPQRSTVWY*')
                if invalid_chars:
                    st.warning(f"Invalid amino acid codes found: {', '.join(invalid_chars)}. These will be skipped.")
                    clean_protein = ''.join(c for c in clean_protein if c in 'ACDEFGHIKLMNPQRSTVWY*')
                
                if clean_protein:
                    with st.spinner("Reverse translating sequence..."):
                        if organism:
                            dna_result = reverse_translate_to_dna(clean_protein, organism, reverse_method)
                        else:
                            # Simple reverse translation without optimization
                            dna_result = reverse_translate_to_dna(clean_protein, 'E. coli', 'most_frequent')
                        
                        st.success("✅ Reverse translation completed successfully!")
                        
                        # Display results
                        st.markdown('<div class="translation-section-header">🔄 Reverse Translation Result</div>', unsafe_allow_html=True)
                        
                        col_result1, col_result2 = st.columns(2)
                        
                        with col_result1:
                            st.markdown("**Original Protein**")
                            st.markdown(f'<div class="translation-sequence-display">{clean_protein}</div>', unsafe_allow_html=True)
                            
                            copy_protein_html = create_copy_button_html(clean_protein, "Protein", "protein2")
                            st.components.v1.html(copy_protein_html, height=80)
                        
                        with col_result2:
                            st.markdown("**Reverse Translated DNA**")
                            st.markdown(f'<div class="translation-sequence-display">{dna_result}</div>', unsafe_allow_html=True)
                            
                            copy_dna_html = create_copy_button_html(dna_result, "DNA", "dna2")
                            st.components.v1.html(copy_dna_html, height=80)
                        
                        # Analysis metrics
                        st.markdown('<div class="translation-section-header">📊 Reverse Translation Metrics</div>', unsafe_allow_html=True)
                        
                        st.markdown('<div class="translation-results-card">', unsafe_allow_html=True)
                        
                        col_metric1, col_metric2, col_metric3, col_metric4 = st.columns(4)
                        
                        with col_metric1:
                            st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                            st.metric("Protein Length", f"{len(clean_protein)} aa")
                            st.markdown('</div>', unsafe_allow_html=True)
                        
                        with col_metric2:
                            st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                            st.metric("DNA Length", f"{len(dna_result)} bp")
                            st.markdown('</div>', unsafe_allow_html=True)
                        
                        with col_metric3:
                            st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                            st.metric("GC Content", f"{calculate_gc_content(dna_result):.1f}%")
                            st.markdown('</div>', unsafe_allow_html=True)
                        
                        with col_metric4:
                            st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                            st.metric("Organism", organism if organism else "None")
                            st.markdown('</div>', unsafe_allow_html=True)
                        
                        st.markdown('</div>', unsafe_allow_html=True)
                        
                        # Verify translation
                        verify_protein = translate_sequence(dna_result, 0, False)
                        
                        if verify_protein == clean_protein:
                            st.markdown('<div class="translation-info-box">✅ Verification passed: Reverse translated DNA produces the original protein!</div>', unsafe_allow_html=True)
                        else:
                            st.markdown('<div class="translation-warning-box">⚠️ Verification note: Some differences may occur due to stop codon handling.</div>', unsafe_allow_html=True)
                        
                        # Download options
                        st.markdown('<div class="translation-section-header">💾 Export Options</div>', unsafe_allow_html=True)
                        
                        col_download1, col_download2 = st.columns(2)
                        
                        with col_download1:
                            protein_fasta = f">Protein_sequence\n{clean_protein}"
                            st.download_button(
                                "📄 Download Protein FASTA",
                                protein_fasta,
                                "protein_sequence.fasta",
                                "text/plain"
                            )
                        
                        with col_download2:
                            organism_name = organism.replace(' ', '_') if organism else 'no_optimization'
                            dna_fasta = f">Reverse_translated_DNA_{organism_name}\n{dna_result}"
                            st.download_button(
                                "📄 Download DNA FASTA",
                                dna_fasta,
                                f"reverse_translated_{organism_name}.fasta",
                                "text/plain",
                                type="primary"
                            )
    
    with tab3:
        st.markdown('<div class="translation-section-header">🔤 Amino Acid Code Conversion</div>', unsafe_allow_html=True)
        
        # Three-letter to one-letter conversion
        st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
        st.markdown("**Three-letter → One-letter Conversion**")
        
        three_letter_input = st.text_input(
            "Enter three-letter amino acid codes (space-separated)",
            placeholder="MET LYS TRP VAL THR PHE ILE SER LEU LEU LEU",
            help="Enter amino acid codes separated by spaces (e.g., MET LYS TRP)"
        )
        
        if st.button("Convert to One-letter", key="three_to_one"):
            if three_letter_input.strip():
                one_letter_result = convert_three_to_one(three_letter_input)
                
                col_conv1, col_conv2 = st.columns(2)
                
                with col_conv1:
                    st.markdown("**Three-letter Input**")
                    st.markdown(f'<div class="translation-sequence-display">{three_letter_input}</div>', unsafe_allow_html=True)
                
                with col_conv2:
                    st.markdown("**One-letter Result**")
                    st.markdown(f'<div class="translation-sequence-display">{one_letter_result}</div>', unsafe_allow_html=True)
                    
                    copy_one_html = create_copy_button_html(one_letter_result, "One-letter", "one1")
                    st.components.v1.html(copy_one_html, height=80)
            else:
                st.error("Please enter three-letter amino acid codes.")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # One-letter to three-letter conversion
        st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
        st.markdown("**One-letter → Three-letter Conversion**")
        
        one_letter_input = st.text_input(
            "Enter one-letter amino acid codes",
            placeholder="MKWVTFISLLL",
            help="Enter single-letter amino acid codes (e.g., MKWVTFISLLL)"
        )
        
        if st.button("Convert to Three-letter", key="one_to_three"):
            if one_letter_input.strip():
                three_letter_result = convert_to_three_letter(one_letter_input)
                
                col_conv1, col_conv2 = st.columns(2)
                
                with col_conv1:
                    st.markdown("**One-letter Input**")
                    st.markdown(f'<div class="translation-sequence-display">{one_letter_input}</div>', unsafe_allow_html=True)
                
                with col_conv2:
                    st.markdown("**Three-letter Result**")
                    st.markdown(f'<div class="translation-sequence-display">{three_letter_result}</div>', unsafe_allow_html=True)
                    
                    copy_three_html = create_copy_button_html(three_letter_result, "Three-letter", "three1")
                    st.components.v1.html(copy_three_html, height=80)
            else:
                st.error("Please enter one-letter amino acid codes.")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Amino acid reference table
        st.markdown('<div class="translation-section-header">📚 Amino Acid Reference Table</div>', unsafe_allow_html=True)
        
        aa_ref_data = []
        for one, props in AA_PROPERTIES.items():
            if one != '*':
                three = ONE_TO_THREE.get(one, 'XXX')
                aa_ref_data.append({
                    'One-letter': one,
                    'Three-letter': three,
                    'Name': props['name'],
                    'MW (Da)': f"{props['mw']:.1f}",
                    'Hydrophobic': '✓' if props['hydrophobic'] else '✗',
                    'Charged': '✓' if props['charged'] else '✗'
                })
        
        aa_ref_df = pd.DataFrame(aa_ref_data)
        st.dataframe(aa_ref_df, use_container_width=True)
    
    with tab4:
        st.markdown('<div class="translation-section-header">📊 Sequence Analysis Tools</div>', unsafe_allow_html=True)
        
        analysis_type = st.selectbox(
            "**Choose Analysis Type**",
            ["Protein Properties", "DNA Analysis", "Codon Usage Comparison"]
        )
        
        if analysis_type == "Protein Properties":
            st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
            st.markdown("**Protein Sequence Analysis**")
            
            protein_for_analysis = st.text_area(
                "Enter protein sequence for analysis",
                height=100,
                placeholder="MKWVTFISLLLLFSSAYSRGVFRRD..."
            )
            st.markdown('</div>', unsafe_allow_html=True)
            
            if protein_for_analysis.strip():
                properties = analyze_protein_properties(protein_for_analysis.strip().upper())
                
                # Display metrics
                st.markdown('<div class="translation-results-card">', unsafe_allow_html=True)
                
                col_prop1, col_prop2, col_prop3, col_prop4 = st.columns(4)
                
                with col_prop1:
                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                    st.metric("Length", f"{properties['length']} aa")
                    st.markdown('</div>', unsafe_allow_html=True)
                
                with col_prop2:
                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                    st.metric("Molecular Weight", f"{properties['molecular_weight']:.1f} Da")
                    st.markdown('</div>', unsafe_allow_html=True)
                
                with col_prop3:
                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                    st.metric("Hydrophobic %", f"{properties.get('hydrophobic_percent', 0):.1f}%")
                    st.markdown('</div>', unsafe_allow_html=True)
                
                with col_prop4:
                    st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                    st.metric("Charged %", f"{properties.get('charged_percent', 0):.1f}%")
                    st.markdown('</div>', unsafe_allow_html=True)
                
                st.markdown('</div>', unsafe_allow_html=True)
                
                # Composition analysis
                if len(protein_for_analysis.strip()) > 1:
                    composition_fig = create_amino_acid_composition_chart(protein_for_analysis.strip().upper())
                    st.plotly_chart(composition_fig, use_container_width=True)
                
                # Detailed composition table
                composition_data = []
                for aa, count in properties['composition'].items():
                    if aa in AA_PROPERTIES:
                        props = AA_PROPERTIES[aa]
                        composition_data.append({
                            'Amino Acid': aa,
                            'Name': props['name'],
                            'Count': count,
                            'Percentage': f"{(count/properties['length']*100):.1f}%",
                            'MW Contribution': f"{(count * props['mw']):.1f} Da"
                        })
                
                comp_df = pd.DataFrame(composition_data).sort_values('Count', ascending=False)
                st.dataframe(comp_df, use_container_width=True)
        
        elif analysis_type == "DNA Analysis":
            st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
            st.markdown("**DNA Sequence Analysis**")
            
            dna_for_analysis = st.text_area(
                "Enter DNA sequence for analysis",
                height=100,
                placeholder="ATGAAGTGGGTAACCTTCATC..."
            )
            st.markdown('</div>', unsafe_allow_html=True)
            
            if dna_for_analysis.strip():
                clean_dna_analysis = clean_dna_sequence(dna_for_analysis)
                
                if clean_dna_analysis:
                    # Basic metrics
                    st.markdown('<div class="translation-results-card">', unsafe_allow_html=True)
                    
                    col_dna1, col_dna2, col_dna3, col_dna4 = st.columns(4)
                    
                    with col_dna1:
                        st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                        st.metric("Length", f"{len(clean_dna_analysis)} bp")
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    with col_dna2:
                        st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                        st.metric("GC Content", f"{calculate_gc_content(clean_dna_analysis):.1f}%")
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    with col_dna3:
                        st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                        reverse_comp = reverse_complement(clean_dna_analysis)
                        st.metric("Reverse Complement", "Available")
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    with col_dna4:
                        st.markdown('<div class="translation-metric-container">', unsafe_allow_html=True)
                        codon_count = len(clean_dna_analysis) // 3
                        st.metric("Complete Codons", codon_count)
                        st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                    # Nucleotide composition
                    nt_composition = Counter(clean_dna_analysis)
                    
                    fig = go.Figure()
                    
                    nucleotides = ['A', 'T', 'G', 'C']
                    counts = [nt_composition.get(nt, 0) for nt in nucleotides]
                    colors = ['#ff4b4b', '#ff6b6b', '#ff8a80', '#ffcdd2']
                    
                    fig.add_trace(go.Bar(
                        x=nucleotides,
                        y=counts,
                        marker_color=colors,
                        text=[f"{c}" for c in counts],
                        textposition='auto'
                    ))
                    
                    fig.update_layout(
                        title="Nucleotide Composition",
                        xaxis_title="Nucleotide",
                        yaxis_title="Count",
                        plot_bgcolor='#f8f9ff',
                        paper_bgcolor='white',
                        height=400
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Show reverse complement
                    st.markdown('<div class="translation-section-header">🔄 Reverse Complement</div>', unsafe_allow_html=True)
                    
                    col_rc1, col_rc2 = st.columns(2)
                    
                    with col_rc1:
                        st.markdown("**Original**")
                        st.markdown(f'<div class="translation-sequence-display">{clean_dna_analysis}</div>', unsafe_allow_html=True)
                    
                    with col_rc2:
                        st.markdown("**Reverse Complement**")
                        st.markdown(f'<div class="translation-sequence-display">{reverse_comp}</div>', unsafe_allow_html=True)
                        
                        copy_rc_html = create_copy_button_html(reverse_comp, "Reverse Complement", "rc1")
                        st.components.v1.html(copy_rc_html, height=80)
        
        elif analysis_type == "Codon Usage Comparison":
            st.markdown('<div class="translation-input-card">', unsafe_allow_html=True)
            st.markdown("**Codon Usage Comparison Between Organisms**")
            
            selected_organisms = st.multiselect(
                "Select organisms to compare",
                list(CODON_USAGE_TABLES.keys()),
                default=list(CODON_USAGE_TABLES.keys())[:2]
            )
            st.markdown('</div>', unsafe_allow_html=True)
            
            if len(selected_organisms) >= 2:
                # Create comparison visualization
                comparison_data = []
                for org in selected_organisms:
                    for aa, codons in CODON_USAGE_TABLES[org].items():
                        for codon, frequency in codons.items():
                            comparison_data.append({
                                'Organism': org,
                                'Amino Acid': aa,
                                'Codon': codon,
                                'Frequency': frequency
                            })
                
                comp_df = pd.DataFrame(comparison_data)
                
                # Show comparison for a specific amino acid
                aa_to_compare = st.selectbox(
                    "Select amino acid to compare",
                    sorted(set(comp_df['Amino Acid']))
                )
                
                filtered_df = comp_df[comp_df['Amino Acid'] == aa_to_compare]
                
                fig = px.bar(
                    filtered_df,
                    x='Codon',
                    y='Frequency',
                    color='Organism',
                    title=f"Codon Usage Comparison for {aa_to_compare}",
                    color_discrete_sequence=['#ff4b4b', '#667eea', '#4CAF50']
                )
                
                fig.update_layout(
                    plot_bgcolor='#f8f9ff',
                    paper_bgcolor='white',
                    height=400
                )
                
                st.plotly_chart(fig, use_container_width=True)
                
                # Show detailed table
                pivot_df = filtered_df.pivot(index='Codon', columns='Organism', values='Frequency')
                st.dataframe(pivot_df, use_container_width=True)

if __name__ == "__main__":
    main()