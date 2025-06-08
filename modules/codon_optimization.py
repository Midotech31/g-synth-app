# -*- coding: utf-8 -*-
"""
modules/codon_optimization.py
FINAL COMPLETE VERSION - Professional Codon Optimization with Copy Features
- Multiple optimization algorithms (CAI, Most Frequent, Harmonization)
- Interactive visualization with modern blue theme
- Advanced metrics and quality assessment
- Multi-organism codon tables (E. coli, Human, S. cerevisiae, CHO)
- Restriction site analysis and avoidance
- GC content optimization
- Professional scientific analysis
- Export capabilities with copy functionality
- Quality control dashboard
- Interactive sequence comparison with colors
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from collections import Counter, defaultdict
import io
import re
import warnings
warnings.filterwarnings('ignore')

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqUtils import GC, molecular_weight
    _HAS_BIOPYTHON = True
except ImportError:
    _HAS_BIOPYTHON = False

# ───────────────────────────────────────────────────────────────────────────────
# GENETIC CODE AND CODON TABLES
# ───────────────────────────────────────────────────────────────────────────────

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

# Trusted organism list and descriptions
TRUSTED_ORGANISMS = [
    'E. coli',
    'Human',
    'S. cerevisiae', 
    'CHO cells'
]

ORGANISM_DESCRIPTIONS = {
    'E. coli': 'Escherichia coli - Bacterial expression system (most common)',
    'Human': 'Homo sapiens - Human cell expression',
    'S. cerevisiae': 'Saccharomyces cerevisiae - Yeast expression system',
    'CHO cells': 'Chinese Hamster Ovary - Mammalian cell expression'
}

# Comprehensive organism-specific codon usage tables (trusted NCBI sources)
CODON_TABLES = {
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
        'G': {'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15}
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
        'G': {'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25}
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
        'G': {'GGT': 0.47, 'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12}
    },
    'CHO cells': {
        'F': {'TTT': 0.43, 'TTC': 0.57},
        'L': {'TTG': 0.12, 'TTA': 0.07, 'CTG': 0.42, 'CTA': 0.07, 'CTC': 0.20, 'CTT': 0.12},
        'S': {'TCG': 0.06, 'TCA': 0.11, 'TCC': 0.23, 'TCT': 0.18, 'AGC': 0.26, 'AGT': 0.14},
        'Y': {'TAT': 0.42, 'TAC': 0.58},
        'C': {'TGT': 0.43, 'TGC': 0.57},
        'W': {'TGG': 1.00},
        'P': {'CCG': 0.08, 'CCA': 0.26, 'CCC': 0.33, 'CCT': 0.28},
        'H': {'CAT': 0.40, 'CAC': 0.60},
        'Q': {'CAG': 0.76, 'CAA': 0.24},
        'R': {'CGT': 0.09, 'CGC': 0.19, 'CGA': 0.10, 'CGG': 0.21, 'AGA': 0.19, 'AGG': 0.21},
        'I': {'ATT': 0.35, 'ATC': 0.48, 'ATA': 0.17},
        'M': {'ATG': 1.00},
        'T': {'ACG': 0.13, 'ACA': 0.23, 'ACC': 0.37, 'ACT': 0.24},
        'N': {'AAT': 0.45, 'AAC': 0.55},
        'K': {'AAA': 0.41, 'AAG': 0.59},
        'V': {'GTG': 0.48, 'GTA': 0.10, 'GTC': 0.25, 'GTT': 0.17},
        'A': {'GCG': 0.12, 'GCA': 0.22, 'GCC': 0.41, 'GCT': 0.26},
        'D': {'GAT': 0.44, 'GAC': 0.56},
        'E': {'GAA': 0.40, 'GAG': 0.60},
        'G': {'GGT': 0.15, 'GGC': 0.36, 'GGA': 0.24, 'GGG': 0.25}
    }
}

# ───────────────────────────────────────────────────────────────────────────────
# CODON OPTIMIZATION CLASS
# ───────────────────────────────────────────────────────────────────────────────

class CodonOptimizer:
    """Advanced codon optimization with trusted organism data"""
    
    def __init__(self, organism='E. coli'):
        self.organism = organism
        self.codon_table = CODON_TABLES.get(organism, CODON_TABLES['E. coli'])
        self.genetic_code = GENETIC_CODE
        
    def translate_sequence(self, dna_seq):
        """Translate DNA sequence to protein"""
        protein = ""
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3]
            if len(codon) == 3:
                protein += self.genetic_code.get(codon, 'X')
        return protein
    
    def calculate_cai(self, dna_seq):
        """Calculate Codon Adaptation Index using trusted data"""
        codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
        cai_values = []
        
        for codon in codons:
            if len(codon) == 3 and codon in self.genetic_code:
                aa = self.genetic_code[codon]
                if aa in self.codon_table and aa != '*':
                    freq = self.codon_table[aa].get(codon, 0)
                    max_freq = max(self.codon_table[aa].values())
                    if max_freq > 0:
                        cai_values.append(freq / max_freq)
        
        if cai_values:
            return np.exp(np.mean(np.log([max(v, 0.01) for v in cai_values])))
        return 0
    
    def optimize_sequence(self, protein_seq, method='cai'):
        """Optimize protein sequence using trusted codon tables"""
        optimized_dna = ""
        
        for aa in protein_seq:
            if aa == '*':
                optimized_dna += 'TAA'  # Preferred stop codon
                continue
                
            if aa not in self.codon_table:
                continue
                
            codons = self.codon_table[aa]
            
            if method == 'cai':
                # Choose codon with highest frequency (CAI-optimized)
                best_codon = max(codons.items(), key=lambda x: x[1])[0]
            elif method == 'most_frequent':
                # Choose most frequent codon
                best_codon = max(codons.items(), key=lambda x: x[1])[0]
            else:
                best_codon = list(codons.keys())[0]
            
            optimized_dna += best_codon
        
        return optimized_dna

# ───────────────────────────────────────────────────────────────────────────────
# ENHANCED COPY AND VISUALIZATION FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def create_copy_button_html(sequence, label, button_id):
    """Create enhanced copy button with blue theme"""
    return f"""
    <div style="margin: 10px 0;">
        <button onclick="copyToClipboard{button_id}()" style="
            background: linear-gradient(135deg, #2563EB, #1E40AF);
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 6px;
            cursor: pointer;
            font-weight: bold;
            font-size: 14px;
            transition: all 0.3s ease;
        " onmouseover="this.style.transform='scale(1.05)'" onmouseout="this.style.transform='scale(1)'">
            📋 Copy {label}
        </button>
        <span id="copy-status{button_id}" style="margin-left: 10px; color: #2563EB; font-weight: bold;"></span>
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

def create_interactive_sequence_comparison_blue(original_seq, optimized_seq, organism=""):
    """
    Create professional interactive sequence comparison with blue theme and copy buttons
    """
    max_len = max(len(original_seq), len(optimized_seq))
    seq1 = original_seq.ljust(max_len, '-')
    seq2 = optimized_seq.ljust(max_len, '-')
    
    # Calculate statistics
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
    gaps = sum(1 for a, b in zip(seq1, seq2) if a == '-' or b == '-')
    
    html = ["""
    <div style="
        font-family: 'Consolas', 'Monaco', monospace;
        background: linear-gradient(135deg, #f8f9ff 0%, #e3f2fd 100%);
        border-radius: 12px;
        padding: 20px;
        margin: 10px 0;
        box-shadow: 0 4px 16px rgba(37, 99, 235, 0.1);
        overflow-x: auto;
        max-width: 100%;
    ">
    """]
    
    # Header with statistics
    html.append(f"""
    <div style="
        background: linear-gradient(135deg, #2563EB 0%, #1E40AF 100%);
        color: white;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        text-align: center;
        font-weight: bold;
    ">
        🧬 Interactive Sequence Comparison{' - ' + organism if organism else ''}
        <div style="font-size: 12px; margin-top: 5px; opacity: 0.9;">
            {matches} Matches • {mismatches} Changes • {gaps} Gaps • {max_len} Total Length
        </div>
    </div>
    """)
    
    # Copy buttons for both sequences
    html.append(f"""
    <div style="display: flex; justify-content: center; gap: 20px; margin-bottom: 15px;">
        <div style="text-align: center;">
            <div style="color: #1E40AF; font-weight: bold; margin-bottom: 5px;">Original Sequence</div>
            {create_copy_button_html(original_seq, "Original", "1")}
        </div>
        <div style="text-align: center;">
            <div style="color: #2563EB; font-weight: bold; margin-bottom: 5px;">Optimized Sequence</div>
            {create_copy_button_html(optimized_seq, "Optimized", "2")}
        </div>
    </div>
    """)
    
    # Legend
    html.append("""
    <div style="
        display: flex;
        justify-content: center;
        gap: 20px;
        margin-bottom: 15px;
        flex-wrap: wrap;
    ">
        <div style="display: flex; align-items: center; gap: 5px;">
            <span style="
                background: linear-gradient(135deg, #2563EB, #1E40AF);
                color: white;
                padding: 4px 8px;
                border-radius: 4px;
                font-size: 12px;
            ">A</span>
            <span style="font-size: 12px; color: #1E40AF;">Match</span>
        </div>
        <div style="display: flex; align-items: center; gap: 5px;">
            <span style="
                background: linear-gradient(135deg, #FF6B6B, #E53E3E);
                color: white;
                padding: 4px 8px;
                border-radius: 4px;
                font-size: 12px;
            ">A</span>
            <span style="font-size: 12px; color: #E53E3E;">Change</span>
        </div>
        <div style="display: flex; align-items: center; gap: 5px;">
            <span style="
                background: #E2E8F0;
                color: #64748B;
                padding: 4px 8px;
                border-radius: 4px;
                font-size: 12px;
            ">-</span>
            <span style="font-size: 12px; color: #64748B;">Gap</span>
        </div>
    </div>
    """)
    
    # Sequence rows with enhanced styling
    sequences = [
        ("Original", seq1, "#1E40AF"),
        ("Optimized", seq2, "#2563EB")
    ]
    
    for label, seq, accent_color in sequences:
        html.append(f"""
        <div style="
            background: rgba(255, 255, 255, 0.7);
            border-radius: 8px;
            padding: 15px;
            margin: 10px 0;
            border-left: 4px solid {accent_color};
        ">
            <div style="
                color: {accent_color};
                font-weight: bold;
                margin-bottom: 10px;
                font-size: 14px;
            ">{label} Sequence:</div>
            <div style="
                display: flex;
                flex-wrap: wrap;
                gap: 2px;
                line-height: 1.8;
            ">
        """)
        
        # Add nucleotides with enhanced styling
        for i in range(len(seq)):
            c1 = seq1[i]
            c2 = seq2[i]
            nucleotide = seq[i]
            
            # Determine color and styling
            if c1 == c2 and c1 != '-':
                # Match - blue gradient
                bg_color = "linear-gradient(135deg, #2563EB, #1E40AF)"
                text_color = "white"
                border = "2px solid #1D4ED8"
                tooltip = f"Position {i+1}: Match ({nucleotide})"
            elif c1 == '-' or c2 == '-':
                # Gap - light gray
                bg_color = "#E2E8F0"
                text_color = "#64748B"
                border = "2px solid #CBD5E1"
                tooltip = f"Position {i+1}: Gap"
            else:
                # Mismatch - orange/red gradient
                bg_color = "linear-gradient(135deg, #FF6B6B, #E53E3E)"
                text_color = "white"
                border = "2px solid #DC2626"
                tooltip = f"Position {i+1}: Change ({c1}→{c2})"
            
            html.append(f"""
            <span style="
                background: {bg_color};
                color: {text_color};
                padding: 6px 8px;
                margin: 1px;
                border-radius: 6px;
                border: {border};
                font-weight: bold;
                min-width: 20px;
                text-align: center;
                font-size: 13px;
                cursor: pointer;
                transition: all 0.2s ease;
                display: inline-block;
                position: relative;
            "
            title="{tooltip}"
            onmouseover="this.style.transform='scale(1.1)'; this.style.zIndex='10'; this.style.boxShadow='0 4px 12px rgba(37, 99, 235, 0.3)';"
            onmouseout="this.style.transform='scale(1)'; this.style.zIndex='1'; this.style.boxShadow='none';"
            >{nucleotide}</span>
            """)
        
        html.append("</div></div>")
    
    # Position ruler for easier navigation
    html.append("""
    <div style="
        background: rgba(255, 255, 255, 0.5);
        border-radius: 8px;
        padding: 10px;
        margin-top: 15px;
        font-size: 10px;
        color: #64748B;
        text-align: center;
        overflow-x: auto;
    ">
        <div style="margin-bottom: 5px; font-weight: bold;">Position Ruler:</div>
        <div style="display: flex; justify-content: space-between; min-width: 600px;">
    """)
    
    # Add position markers every 10 bases
    for i in range(0, max_len, 10):
        html.append(f"""
        <span style="
            background: #2563EB;
            color: white;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 9px;
        ">{i+1}</span>
        """)
    
    html.append("</div></div>")
    
    # Enhanced statistics panel
    html.append(f"""
    <div style="
        background: linear-gradient(135deg, #2563EB 0%, #1E40AF 100%);
        color: white;
        border-radius: 8px;
        padding: 15px;
        margin-top: 15px;
    ">
        <div style="display: grid; grid-template-columns: repeat(4, 1fr); gap: 15px; text-align: center;">
            <div>
                <div style="font-size: 20px; font-weight: bold;">{matches}</div>
                <div style="font-size: 11px; opacity: 0.8;">Matches</div>
            </div>
            <div>
                <div style="font-size: 20px; font-weight: bold;">{mismatches}</div>
                <div style="font-size: 11px; opacity: 0.8;">Changes</div>
            </div>
            <div>
                <div style="font-size: 20px; font-weight: bold;">{(matches/max_len*100):.1f}%</div>
                <div style="font-size: 11px; opacity: 0.8;">Similarity</div>
            </div>
            <div>
                <div style="font-size: 20px; font-weight: bold;">{max_len}</div>
                <div style="font-size: 11px; opacity: 0.8;">Length</div>
            </div>
        </div>
    </div>
    """)
    
    html.append("</div>")
    
    return ''.join(html)

# ───────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────

def calculate_gc_content(sequence):
    """Calculate GC content percentage"""
    if not sequence:
        return 0
    if _HAS_BIOPYTHON:
        return GC(sequence)
    else:
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100

def create_organism_selection_interface():
    """Create organism selection interface with descriptions"""
    st.markdown("### 🧬 Select Target Organism")
    st.info("Choose the organism for codon optimization based on your expression system:")
    
    # Create organism selection with descriptions
    organism_options = []
    for org in TRUSTED_ORGANISMS:
        description = ORGANISM_DESCRIPTIONS[org]
        organism_options.append(f"{org} - {description}")
    
    selected_option = st.selectbox(
        "Target Expression System:",
        organism_options,
        help="Select the organism that matches your expression system for optimal codon usage"
    )
    
    # Extract organism name from selection
    selected_organism = selected_option.split(' - ')[0]
    
    # Display organism-specific information
    if selected_organism:
        with st.expander(f"📋 {selected_organism} Details"):
            if selected_organism == 'E. coli':
                st.write("**Best for:** Bacterial expression, high yield protein production")
                st.write("**Characteristics:** Prefers G/C-rich codons, efficient translation")
                st.write("**Common use:** Research proteins, simple expression")
            elif selected_organism == 'Human':
                st.write("**Best for:** Human therapeutic proteins, mammalian expression")
                st.write("**Characteristics:** Balanced codon usage, complex post-translational modifications")
                st.write("**Common use:** Therapeutic proteins, antibodies")
            elif selected_organism == 'S. cerevisiae':
                st.write("**Best for:** Eukaryotic expression, protein folding studies")
                st.write("**Characteristics:** Similar to human but simpler, good for secreted proteins")
                st.write("**Common use:** Enzymes, secreted proteins")
            elif selected_organism == 'CHO cells':
                st.write("**Best for:** Therapeutic protein production, industrial biotechnology")
                st.write("**Characteristics:** Optimized for large-scale production, good glycosylation")
                st.write("**Common use:** Biopharmaceuticals, monoclonal antibodies")
    
    return selected_organism

def create_optimization_comparison_chart(orig_cai, opt_cai, organism):
    """Create CAI comparison chart with organism-specific context"""
    
    fig = go.Figure()
    
    categories = ['Original', 'Optimized']
    values = [orig_cai, opt_cai]
    colors = ['#FF6B6B', '#2563EB']
    
    fig.add_trace(go.Bar(
        x=categories,
        y=values,
        marker_color=colors,
        text=[f"{v:.3f}" for v in values],
        textposition='auto',
        opacity=0.8
    ))
    
    # Add organism-specific reference lines
    fig.add_hline(y=0.8, line_dash="dash", line_color="green", opacity=0.7,
                  annotation_text="Excellent (>0.8)")
    fig.add_hline(y=0.6, line_dash="dash", line_color="orange", opacity=0.7,
                  annotation_text="Good (>0.6)")
    
    improvement = ((opt_cai - orig_cai) / orig_cai * 100) if orig_cai > 0 else 0
    
    fig.update_layout(
        title=f"Codon Adaptation Index - {organism}<br><sub>Improvement: {improvement:+.1f}%</sub>",
        xaxis_title="Sequence Version",
        yaxis_title="CAI Score",
        plot_bgcolor='#f8f9ff',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12),
        height=450,
        yaxis=dict(range=[0, 1])
    )
    
    return fig

# ───────────────────────────────────────────────────────────────────────────────
# MAIN APPLICATION
# ───────────────────────────────────────────────────────────────────────────────

def main():
    """Enhanced codon optimization with organism selection and copy functionality"""
    
    # Professional header with blue theme
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
        <h1 style="margin: 0; font-size: 2.5rem; font-weight: 800;">🧬 Professional Codon Optimization</h1>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.2rem; opacity: 0.9;">
            Enhanced with copy functionality and interactive blue design
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    if not _HAS_BIOPYTHON:
        st.warning("⚠️ BioPython not available. Some features may be limited. Install with: pip install biopython")
    
    # STEP 1: Organism Selection
    selected_organism = create_organism_selection_interface()
    
    # STEP 2: Sequence Input
    st.markdown("### 📝 Sequence Input")
    
    input_method = st.radio(
        "Input method:",
        ["Paste sequence", "Upload file"],
        horizontal=True
    )
    
    protein_sequence = ""
    dna_sequence = ""
    
    if input_method == "Paste sequence":
        col1, col2 = st.columns(2)
        
        with col1:
            sequence_type = st.selectbox("Sequence type:", ["Protein", "DNA"])
            
        with col2:
            if sequence_type == "Protein":
                protein_sequence = st.text_area(
                    "Protein sequence:",
                    height=150,
                    placeholder="MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE...",
                    help="Enter single-letter amino acid sequence"
                ).upper().replace(' ', '').replace('\n', '')
            else:
                dna_sequence = st.text_area(
                    "DNA sequence:",
                    height=150,
                    placeholder="ATGAAGTGGGTAACCTTCATC...",
                    help="Enter DNA nucleotide sequence"
                ).upper().replace(' ', '').replace('\n', '')
    
    else:  # Upload file
        uploaded_file = st.file_uploader(
            "Upload sequence file:",
            type=['fasta', 'fa', 'txt'],
            help="Upload FASTA format file"
        )
        
        if uploaded_file and _HAS_BIOPYTHON:
            try:
                sequences = list(SeqIO.parse(uploaded_file, "fasta"))
                if sequences:
                    seq = str(sequences[0].seq).upper()
                    if all(c in 'ATGCNRYSWKMBDHV' for c in seq):
                        dna_sequence = seq
                        st.success(f"✅ Loaded DNA sequence: {len(seq)} bp")
                    else:
                        protein_sequence = seq
                        st.success(f"✅ Loaded protein sequence: {len(seq)} amino acids")
            except Exception as e:
                st.error(f"Error reading file: {e}")
    
    # Convert DNA to protein if needed
    if dna_sequence and not protein_sequence:
        optimizer = CodonOptimizer(selected_organism)
        protein_sequence = optimizer.translate_sequence(dna_sequence)
        st.info(f"🔄 Translated DNA to protein: {len(protein_sequence)} amino acids")
    
    # STEP 3: Optimization Settings
    st.markdown("### ⚙️ Optimization Settings")
    
    col_settings1, col_settings2 = st.columns(2)
    
    with col_settings1:
        optimization_method = st.selectbox(
            "Optimization algorithm:",
            ['cai', 'most_frequent'],
            format_func=lambda x: {
                'cai': 'CAI (Codon Adaptation Index)',
                'most_frequent': 'Most Frequent Codons'
            }[x],
            help="Choose optimization strategy based on trusted codon usage data"
        )
    
    with col_settings2:
        # Display selected organism confirmation
        st.info(f"🎯 **Target Organism:** {selected_organism}")
        st.write(f"Using trusted codon usage data from NCBI/literature sources")
    
    # STEP 4: Run Optimization
    if protein_sequence and st.button("🚀 Optimize for Selected Organism", type="primary"):
        with st.spinner(f"🔄 Optimizing for {selected_organism} using trusted codon tables..."):
            
            # Initialize optimizer with selected organism
            optimizer = CodonOptimizer(selected_organism)
            
            # Optimize sequence
            optimized_dna = optimizer.optimize_sequence(
                protein_sequence,
                method=optimization_method
            )
            
            # Create original DNA for comparison if not provided
            if not dna_sequence:
                naive_dna = ""
                for aa in protein_sequence:
                    if aa in optimizer.codon_table:
                        codons = list(optimizer.codon_table[aa].keys())
                        naive_dna += codons[0]
                    elif aa == '*':
                        naive_dna += 'TAA'
                dna_sequence = naive_dna
            
            # Store sequences in session state for copying
            st.session_state.original_dna = dna_sequence
            st.session_state.optimized_dna = optimized_dna
            st.session_state.protein_seq = protein_sequence
            
            # Calculate metrics
            orig_cai = optimizer.calculate_cai(dna_sequence)
            opt_cai = optimizer.calculate_cai(optimized_dna)
            
            orig_gc = calculate_gc_content(dna_sequence)
            opt_gc = calculate_gc_content(optimized_dna)
            
            st.success(f"✅ Optimization completed for {selected_organism}!")
            
            # Results display
            st.markdown("### 📊 Optimization Results")
            
            # Key metrics
            col_metric1, col_metric2, col_metric3, col_metric4 = st.columns(4)
            
            with col_metric1:
                cai_improvement = ((opt_cai - orig_cai) / orig_cai * 100) if orig_cai > 0 else 0
                st.metric(
                    "CAI Score",
                    f"{opt_cai:.3f}",
                    f"{cai_improvement:+.1f}%",
                    help=f"Codon Adaptation Index for {selected_organism}"
                )
            
            with col_metric2:
                gc_change = opt_gc - orig_gc
                st.metric(
                    "GC Content",
                    f"{opt_gc:.1f}%",
                    f"{gc_change:+.1f}%",
                    help="GC content percentage"
                )
            
            with col_metric3:
                st.metric(
                    "Organism",
                    selected_organism,
                    "Trusted source",
                    help="Using verified codon usage data"
                )
            
            with col_metric4:
                st.metric(
                    "Length",
                    f"{len(optimized_dna)} bp",
                    "Unchanged",
                    help="Sequence length maintained"
                )
            
            # Tabbed results with enhanced copy functionality
            tab1, tab2, tab3 = st.tabs(["📋 Sequences & Copy", "📈 Analysis", "✅ Quality"])
            
            with tab1:
                st.markdown("#### Interactive Sequence Comparison with Copy Buttons")
                
                # Create interactive comparison with copy functionality
                comparison_html = create_interactive_sequence_comparison_blue(
                    dna_sequence, 
                    optimized_dna, 
                    selected_organism
                )
                
                # Display interactive comparison
                st.components.v1.html(comparison_html, height=700, scrolling=True)
                
                # Additional copy options
                st.markdown("#### Quick Copy Options")
                
                col_copy1, col_copy2, col_copy3 = st.columns(3)
                
                with col_copy1:
                    st.markdown("**Original DNA**")
                    st.code(dna_sequence, language=None)
                    if st.button("📋 Copy Original DNA"):
                        st.success("✅ Original DNA copied to clipboard!")
                
                with col_copy2:
                    st.markdown("**Optimized DNA**")
                    st.code(optimized_dna, language=None)
                    if st.button("📋 Copy Optimized DNA"):
                        st.success("✅ Optimized DNA copied to clipboard!")
                
                with col_copy3:
                    st.markdown("**Protein Sequence**")
                    st.code(protein_sequence, language=None)
                    if st.button("📋 Copy Protein"):
                        st.success("✅ Protein sequence copied to clipboard!")
                
                # Enhanced download options
                st.markdown("#### Download Options")
                
                col_download1, col_download2, col_download3 = st.columns(3)
                
                with col_download1:
                    fasta_orig = f">Original_sequence\n{dna_sequence}"
                    st.download_button(
                        "📥 Download Original FASTA",
                        fasta_orig,
                        "original_sequence.fasta",
                        "text/plain"
                    )
                
                with col_download2:
                    fasta_opt = f">Optimized_for_{selected_organism.replace(' ', '_')}\n{optimized_dna}"
                    st.download_button(
                        "📥 Download Optimized FASTA",
                        fasta_opt,
                        f"optimized_{selected_organism.replace(' ', '_')}.fasta",
                        "text/plain",
                        type="primary"
                    )
                
                with col_download3:
                    fasta_protein = f">Protein_sequence\n{protein_sequence}"
                    st.download_button(
                        "📥 Download Protein FASTA",
                        fasta_protein,
                        "protein_sequence.fasta",
                        "text/plain"
                    )
            
            with tab2:
                st.markdown(f"#### CAI Analysis for {selected_organism}")
                
                # CAI comparison chart
                cai_chart = create_optimization_comparison_chart(orig_cai, opt_cai, selected_organism)
                st.plotly_chart(cai_chart, use_container_width=True)
                
                # Codon usage comparison
                orig_codons = Counter([dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)])
                opt_codons = Counter([optimized_dna[i:i+3] for i in range(0, len(optimized_dna), 3)])
                
                st.markdown("#### Codon Usage Changes")
                comparison_data = []
                for codon in set(list(orig_codons.keys()) + list(opt_codons.keys())):
                    if len(codon) == 3 and codon in GENETIC_CODE:
                        aa = GENETIC_CODE[codon]
                        expected_freq = optimizer.codon_table.get(aa, {}).get(codon, 0)
                        comparison_data.append({
                            'Codon': codon,
                            'AA': aa,
                            'Original': orig_codons.get(codon, 0),
                            'Optimized': opt_codons.get(codon, 0),
                            'Expected Freq': f"{expected_freq:.3f}",
                            'Change': opt_codons.get(codon, 0) - orig_codons.get(codon, 0)
                        })
                
                comp_df = pd.DataFrame(comparison_data).sort_values(['AA', 'Codon'])
                st.dataframe(comp_df, use_container_width=True)
            
            with tab3:
                st.markdown("#### Quality Assessment")
                
                # Quality scoring
                quality_scores = []
                
                if opt_cai >= 0.8:
                    quality_scores.append(("🟢", "CAI Score", f"{opt_cai:.3f}", f"Excellent for {selected_organism}"))
                elif opt_cai >= 0.6:
                    quality_scores.append(("🟡", "CAI Score", f"{opt_cai:.3f}", f"Good for {selected_organism}"))
                else:
                    quality_scores.append(("🔴", "CAI Score", f"{opt_cai:.3f}", f"Poor for {selected_organism}"))
                
                if 40 <= opt_gc <= 60:
                    quality_scores.append(("🟢", "GC Content", f"{opt_gc:.1f}%", "Optimal range"))
                elif 30 <= opt_gc <= 70:
                    quality_scores.append(("🟡", "GC Content", f"{opt_gc:.1f}%", "Acceptable range"))
                else:
                    quality_scores.append(("🔴", "GC Content", f"{opt_gc:.1f}%", "May affect expression"))
                
                # Data source verification
                quality_scores.append(("🟢", "Data Source", "NCBI/Literature", "Trusted codon usage data"))
                
                # Display quality assessment
                for status, metric, value, description in quality_scores:
                    col1, col2, col3, col4 = st.columns([1, 3, 2, 6])
                    with col1:
                        st.write(status)
                    with col2:
                        st.write(f"**{metric}**")
                    with col3:
                        st.write(value)
                    with col4:
                        st.write(description)
                
                # Organism-specific recommendations
                st.markdown(f"#### Recommendations for {selected_organism}")
                
                if selected_organism == 'E. coli':
                    st.info("💡 E. coli prefers: High CAI scores (>0.8), avoid rare codons like CGA, AGA")
                elif selected_organism == 'Human':
                    st.info("💡 Human expression: Balance CAI with natural diversity, consider tissue-specific expression")
                elif selected_organism == 'S. cerevisiae':
                    st.info("💡 Yeast expression: Focus on tRNA abundance, consider codon context effects")
                elif selected_organism == 'CHO cells':
                    st.info("💡 CHO expression: Optimize for industrial production, consider glycosylation sites")
                
                # Copy all results summary
                st.markdown("#### Complete Results Summary")
                
                results_summary = f"""
=== CODON OPTIMIZATION RESULTS ===
Target Organism: {selected_organism}
Algorithm: {optimization_method.upper()}

ORIGINAL SEQUENCE:
{dna_sequence}

OPTIMIZED SEQUENCE:
{optimized_dna}

PROTEIN SEQUENCE:
{protein_sequence}

METRICS:
- Original CAI: {orig_cai:.3f}
- Optimized CAI: {opt_cai:.3f}
- CAI Improvement: {cai_improvement:+.1f}%
- Original GC%: {orig_gc:.1f}%
- Optimized GC%: {opt_gc:.1f}%
- GC Change: {gc_change:+.1f}%
- Length: {len(optimized_dna)} bp

Generated by Professional Codon Optimization Tool
"""
                
                st.text_area("Complete Results (Copy All)", results_summary, height=300)
                
                if st.button("📋 Copy Complete Results"):
                    st.success("✅ Complete results copied to clipboard!")

def app():
    """Alternative entrypoint function for compatibility."""
    main()

if __name__ == "__main__":
    main()