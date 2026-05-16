# modules/home.py

import streamlit as st

def main():
    # ── 1. Enhanced Logo & Title Layout (side by side, no extra spacing) ──────────
    logo_path = "assets/logo1.png"
    
    # Create columns for logo and title layout
    col_logo, col_title = st.columns([1, 4])
    
    with col_logo:
        try:
            st.image(logo_path, width=120)
        except:
            pass
    
    with col_title:
        st.markdown(
            """
            <div style='padding-top: 10px;'>
                <h1 style='margin: 0; font-size: 2.5rem; font-weight: 700; color: #0b1f3b; line-height: 1.2;'>
                    Welcome to G-Synth 2025.2.0
                </h1>
                <h2 style='margin: 0.3rem 0 0 0; font-size: 1.3rem; font-weight: 500; color: #2563EB; line-height: 1.3;'>
                    The Ultimate Full-Stack Toolkit for Therapeutic Peptide Gene Synthesis
                </h2>
            </div>
            """,
            unsafe_allow_html=True
        )

    # Enhanced description
    st.markdown(
        """
        <div style='margin: 1.5rem 0 2rem 0; max-width: 900px;'>
            <p style='font-size: 1.1rem; color: #374151; line-height: 1.6; margin: 0; text-align: justify;'>
                G-Synth is an advanced software platform that supports the entire workflow of therapeutic peptide gene synthesis — from initial design to final validation. Perform <strong>in silico</strong> gene design, codon optimization, hybridization simulation, cloning validation, reference sequence alignment, and much more — all from a single, streamlined interface.
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

    # ── 2. UPDATED CARDS WITH MORE RELEVANT ICONS ───────────────────
    cards = [
        ("🔬 Small Sequence Design",
         "Reverse complement sequences, calculate GC% and Tₘ, find ORFs, design primers and oligos.",
         "Small Sequence Design"),
        ("🧬 Translation Tools",
         "Convert DNA → protein, reverse-translate protein → DNA, and switch between one/three-letter codes.",
         "Translation / Reverse Transl."),
        ("⚙️ Codon Optimization",
         "Optimize codon usage (E. coli, yeast, CHO, human) – tune GC content and avoid restriction sites.",
         "Codon Optimization"),
        ("🔗 Extended Synthesis",
         "Fragment long DNA into overlapping oligos with sticky ends, verify assembly, export protocols.",
         "Extended Synthesis"),
        ("🧪 Hybridization & Ligation Check",
         "Calculate melting temps for hybridization and check ligation compatibility of fragments.",
         "Hybridization"),
        ("🎯 Primer Generator",
         "Design primers for PCR, sequencing, or site-directed mutagenesis with adjustable parameters.",
         "Primer Generator"),
        ("🔄 Reverse Complement",
         "Instantly get the reverse complement of any DNA sequence.",
         "Reverse Complement"),
        ("✂️ CRISPR sgRNA Designer",
         "Scan for candidate guides, score on-target/off-target (Doench 2016), and link to UCSC genome browser.",
         "CRISPR sgRNA Designer"),
        ("🧫 Plasmid Visualizer & Editor",
         "Upload GenBank/FASTA, annotate features, render maps, and export SVG/PNG/GenBank.",
         "Plasmid Visualizer"),
        ("🧩 Ligation Calculator",
         "Verify overhang complementarity, calculate ΔG/efficiency, and visualize enzyme cut sites.",
         "Ligation Calculator"),
        ("🖥️ In Silico Docking",
         "Perform AI-based docking simulations and visualize predicted binding modes.",
         "In Silico Docking"),
        ("🧠 Functional Prediction",
         "Predict protein function from sequence using machine learning models.",
         "Functional Prediction"),
        ("📊 Alignment Tools",
         "Perform multiple sequence alignments and visualize results.",
         "Alignment Tools"),
        ("⚡ Settings & Preferences",
         "Toggle Dark Mode, manage user settings, and customize your G-Synth workspace.",
         "Settings"),
        ("📖 Help & Guide",
         "Comprehensive user manual, FAQs, and troubleshooting tips.",
         "Help & Guide"),
    ]

    # ── 3. Render cards in rows of 4 for a more compact, classy look ────────────
    cards_per_row = 4
    for row_start in range(0, len(cards), cards_per_row):
        row_cards = cards[row_start:row_start + cards_per_row]
        cols = st.columns(cards_per_row, gap="small")
        for idx, (title, desc, label) in enumerate(row_cards):
            col = cols[idx]
            href = f"/?page={label.replace(' ', '+')}"
            # Enhanced card styling
            card_html = f"""
                <a href="{href}" style="text-decoration:none;">
                  <div class="card"
                       style="
                         border-left: 5px solid #ff4b4b;
                         padding: 0.6rem 0.8rem !important;
                         margin-bottom: 0.75rem;
                         max-width: 280px;
                         cursor: pointer;
                         background: white;
                         border-radius: 8px;
                         box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                         transition: transform 0.2s ease, box-shadow 0.2s ease;
                       "
                       onmouseover="this.style.transform='translateY(-2px)'; this.style.boxShadow='0 4px 12px rgba(0,0,0,0.15)';"
                       onmouseout="this.style.transform='translateY(0px)'; this.style.boxShadow='0 2px 8px rgba(0,0,0,0.1)';">
                      <h3 style="margin-bottom: 0.2rem; color:#0b1f3b; font-size:1rem;">{title}</h3>
                      <p style="font-size:0.85rem; color:#555; margin:0;">{desc}</p>
                  </div>
                </a>
            """
            with col:
                st.markdown(card_html, unsafe_allow_html=True)

    # ── 4. Separator ───────────────────────────────────────────────────────────────
    st.markdown("<hr style='margin-top:1.5rem; margin-bottom:1.5rem; border: 1px solid #e5e7eb;'>", unsafe_allow_html=True)

    # ── 5. ENHANCED Getting Started Section with PROPER VISUAL CARDS ─────────────
    
    # Title and subtitle
    st.markdown("""
    <div style='text-align: center; margin-bottom: 2rem;'>
        <h3 style='color: #0b1f3b; margin: 0 0 0.5rem 0; font-size: 1.8rem; font-weight: 700;'>
            🚀 Getting Started with G-Synth 2025.2.0
        </h3>
        <p style='color: #475569; font-size: 1.1rem; margin: 0; font-weight: 500;'>
            Streamlined workflow for therapeutic peptide gene synthesis
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Create enhanced visual cards using Streamlit columns
    col1, col2 = st.columns(2, gap="large")
    
    with col1:
        # Card 1: In Silico Design
        st.markdown("""
        <div style='
            background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
            padding: 2rem;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            box-shadow: 0 8px 25px rgba(59, 130, 246, 0.1);
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        '>
            <div style='
                position: absolute;
                top: -20px;
                right: -20px;
                width: 60px;
                height: 60px;
                background: linear-gradient(45deg, #3b82f6, #1d4ed8);
                border-radius: 50%;
                opacity: 0.1;
            '></div>
            <div style='position: relative; z-index: 1;'>
                <div style='
                    display: flex;
                    align-items: center;
                    margin-bottom: 1rem;
                    padding-bottom: 0.5rem;
                    border-bottom: 2px solid #e2e8f0;
                '>
                    <span style='font-size: 2rem; margin-right: 0.5rem;'>🧬</span>
                    <h4 style='
                        color: #3b82f6;
                        margin: 0;
                        font-size: 1.25rem;
                        font-weight: 700;
                    '>In Silico Design</h4>
                </div>
                <p style='
                    color: #64748b;
                    font-size: 1rem;
                    margin: 0;
                    line-height: 1.6;
                '>
                    Design peptide-coding DNA sequences with real-time GC%, Tₘ, and ORF analysis.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Card 3: Hybridization Simulation
        st.markdown("""
        <div style='
            background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
            padding: 2rem;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            box-shadow: 0 8px 25px rgba(139, 92, 246, 0.1);
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        '>
            <div style='
                position: absolute;
                top: -20px;
                right: -20px;
                width: 60px;
                height: 60px;
                background: linear-gradient(45deg, #8b5cf6, #7c3aed);
                border-radius: 50%;
                opacity: 0.1;
            '></div>
            <div style='position: relative; z-index: 1;'>
                <div style='
                    display: flex;
                    align-items: center;
                    margin-bottom: 1rem;
                    padding-bottom: 0.5rem;
                    border-bottom: 2px solid #e2e8f0;
                '>
                    <span style='font-size: 2rem; margin-right: 0.5rem;'>🧪</span>
                    <h4 style='
                        color: #8b5cf6;
                        margin: 0;
                        font-size: 1.25rem;
                        font-weight: 700;
                    '>Hybridization Simulation</h4>
                </div>
                <p style='
                    color: #64748b;
                    font-size: 1rem;
                    margin: 0;
                    line-height: 1.6;
                '>
                    Predict oligo binding, secondary structures, and off-target interactions.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Card 5: Reference Sequence Alignment
        st.markdown("""
        <div style='
            background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
            padding: 2rem;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            box-shadow: 0 8px 25px rgba(239, 68, 68, 0.1);
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        '>
            <div style='
                position: absolute;
                top: -20px;
                right: -20px;
                width: 60px;
                height: 60px;
                background: linear-gradient(45deg, #ef4444, #dc2626);
                border-radius: 50%;
                opacity: 0.1;
            '></div>
            <div style='position: relative; z-index: 1;'>
                <div style='
                    display: flex;
                    align-items: center;
                    margin-bottom: 1rem;
                    padding-bottom: 0.5rem;
                    border-bottom: 2px solid #e2e8f0;
                '>
                    <span style='font-size: 2rem; margin-right: 0.5rem;'>📊</span>
                    <h4 style='
                        color: #ef4444;
                        margin: 0;
                        font-size: 1.25rem;
                        font-weight: 700;
                    '>Reference Sequence Alignment</h4>
                </div>
                <p style='
                    color: #64748b;
                    font-size: 1rem;
                    margin: 0;
                    line-height: 1.6;
                '>
                    Align designed sequences against reference genomes and visualize mutations or mismatches.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        # Card 2: Codon Optimization
        st.markdown("""
        <div style='
            background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
            padding: 2rem;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            box-shadow: 0 8px 25px rgba(16, 185, 129, 0.1);
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        '>
            <div style='
                position: absolute;
                top: -20px;
                right: -20px;
                width: 60px;
                height: 60px;
                background: linear-gradient(45deg, #10b981, #059669);
                border-radius: 50%;
                opacity: 0.1;
            '></div>
            <div style='position: relative; z-index: 1;'>
                <div style='
                    display: flex;
                    align-items: center;
                    margin-bottom: 1rem;
                    padding-bottom: 0.5rem;
                    border-bottom: 2px solid #e2e8f0;
                '>
                    <span style='font-size: 2rem; margin-right: 0.5rem;'>⚙️</span>
                    <h4 style='
                        color: #10b981;
                        margin: 0;
                        font-size: 1.25rem;
                        font-weight: 700;
                    '>Codon Optimization</h4>
                </div>
                <p style='
                    color: #64748b;
                    font-size: 1rem;
                    margin: 0;
                    line-height: 1.6;
                '>
                    Optimize sequences for target organisms using built-in codon usage tables and expression tuning.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Card 4: Cloning Validation
        st.markdown("""
        <div style='
            background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
            padding: 2rem;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            box-shadow: 0 8px 25px rgba(245, 158, 11, 0.1);
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        '>
            <div style='
                position: absolute;
                top: -20px;
                right: -20px;
                width: 60px;
                height: 60px;
                background: linear-gradient(45deg, #f59e0b, #d97706);
                border-radius: 50%;
                opacity: 0.1;
            '></div>
            <div style='position: relative; z-index: 1;'>
                <div style='
                    display: flex;
                    align-items: center;
                    margin-bottom: 1rem;
                    padding-bottom: 0.5rem;
                    border-bottom: 2px solid #e2e8f0;
                '>
                    <span style='font-size: 2rem; margin-right: 0.5rem;'>🧫</span>
                    <h4 style='
                        color: #f59e0b;
                        margin: 0;
                        font-size: 1.25rem;
                        font-weight: 700;
                    '>Cloning Validation</h4>
                </div>
                <p style='
                    color: #64748b;
                    font-size: 1rem;
                    margin: 0;
                    line-height: 1.6;
                '>
                    Simulate cloning steps, verify constructs, and preview vector-insert assemblies.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Card 6: All-in-One Interface
        st.markdown("""
        <div style='
            background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
            padding: 2rem;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            box-shadow: 0 8px 25px rgba(6, 182, 212, 0.1);
            margin-bottom: 1.5rem;
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        '>
            <div style='
                position: absolute;
                top: -20px;
                right: -20px;
                width: 60px;
                height: 60px;
                background: linear-gradient(45deg, #06b6d4, #0891b2);
                border-radius: 50%;
                opacity: 0.1;
            '></div>
            <div style='position: relative; z-index: 1;'>
                <div style='
                    display: flex;
                    align-items: center;
                    margin-bottom: 1rem;
                    padding-bottom: 0.5rem;
                    border-bottom: 2px solid #e2e8f0;
                '>
                    <span style='font-size: 2rem; margin-right: 0.5rem;'>🎯</span>
                    <h4 style='
                        color: #06b6d4;
                        margin: 0;
                        font-size: 1.25rem;
                        font-weight: 700;
                    '>All-in-One Interface</h4>
                </div>
                <p style='
                    color: #64748b;
                    font-size: 1rem;
                    margin: 0;
                    line-height: 1.6;
                '>
                    Visualize plasmids, design CRISPR guides, manage settings — all in one seamless workspace.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # ── 6. Footer Ribbon ───────────────────────────────────────────────────────────
    st.markdown(
        """
        <div style="
            background-color: #0b1f3b;
            color: #ffffff;
            padding: 1.5rem 2rem;
            margin-top: 2rem;
            text-align: center;
            border-radius: 4px;
            position: relative;
            overflow: hidden;
        ">
            <div style="
                position: absolute;
                left: -60px;
                top: 50%;
                transform: translateY(-50%) rotate(-45deg);
                width: 200px;
                height: 200px;
                background-color: #ff4b4b;
                z-index: 0;
            "></div>
            <div style="
                position: relative;
                z-index: 1;
                max-width: 800px;
                margin: 0 auto;
            ">
                <h4 style="margin: 0 0 0.5rem 0; font-size:1.2rem; font-weight:600;">👤 Developer</h4>
                <p style="margin: 0.2rem 0 0.5rem 0; font-size:1rem;">
                    <strong>Dr. Mohamed Merzoug</strong> &middot; Genomics Technology Platform, Higher School of Biological Sciences of Oran &middot;
                    <a href="mailto:mohamed.merzoug.essbo@gmail.com" style="color:#ffffff; text-decoration:underline;">mohamed.merzoug.essbo@gmail.com</a>
                </p>
                <p style="margin: 0.5rem 0 0 0; font-size:0.9rem; color:#dddddd;">
                    &copy; 2025 G-Synth Toolkit. All rights reserved.
                </p>
            </div>
        </div>
        """,
        unsafe_allow_html=True
    )

def app():
    """Entry point for modular integration"""
    main()

if __name__ == "__main__":
    main()