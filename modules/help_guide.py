# modules/help_guide.py
# -*- coding: utf-8 -*-

import streamlit as st

def main():
    st.title("G-Synth 2025.2.0 Help & Guide")

    st.markdown(
        """
Welcome to **G-Synth**, your all-in-one gene synthesis toolkit. Use the expandable sections below to explore each module’s features and usage examples.  
If you have questions, expand the relevant section or contact Dr. Mohamed Merzoug at **mohamed.merzoug.essbo@gmail.com**.
"""
    )

    st.markdown("---")

    # 1. Small Sequence Design (SSD)
    with st.expander("1. Small Sequence Design (SSD)", expanded=False):
        st.write(
            """
**Purpose:** Quickly analyze or manipulate short DNA/protein sequences.

#### 1.1 Reverse Complement
- **Function:** Returns the reverse complement of a DNA sequence, preserving ambiguous bases (N, R, Y, etc.).
- **How to use:**
  1. Paste or type a DNA sequence (A/T/C/G/N/R/Y).
  2. Click **Compute Reverse Complement**.
  3. View the cleaned input and its reverse complement.

#### 1.2 GC% Calculator
- **Function:** Computes GC content (percentage) of a DNA sequence.
- **How to use:**
  1. Enter a DNA sequence.
  2. Click **Calculate GC%**.
  3. See total length and the GC percentage.

#### 1.3 Tm (Melting Temperature) Calculator
- **Function:** Uses nearest-neighbor thermodynamic parameters to estimate primer Tₘ.
- **How to use:**
  1. Input primer sequence (DNA).
  2. Specify salt concentration (e.g., 50 mM Na⁺) and oligo concentration (e.g., 0.25 µM).
  3. Click **Compute Tₘ**.
  4. View Tₘ in °C and any warnings for secondary structures or extreme GC%.

#### 1.4 ORF Finder
- **Function:** Scans all six reading frames for open reading frames (ORFs) ≥ a chosen length (default 50 aa).
- **How to use:**
  1. Enter a full DNA sequence (e.g., plasmid or gene).
  2. (Optional) Set minimum ORF length in amino acids.
  3. Click **Find ORFs**.
  4. Inspect a results table listing frame, start/end positions, length, and translated amino acids.

#### 1.5 Basic Codon Optimization
- **Function:** Converts a protein sequence to DNA using a simple codon table.
- **How to use:**
  1. Select **Protein Input** and paste amino acid sequence (one-letter codes).
  2. Choose a target organism (e.g., _E. coli BL21_, _S. cerevisiae_).
  3. Click **Optimize**.
  4. Receive a “quick” optimized DNA (no advanced GC or restriction filtering).
"""
        )

    # 2. Translation & Reverse Translation
    with st.expander("2. Translation & Reverse Translation", expanded=False):
        st.write(
            """
**Purpose:** Convert between DNA and protein with organism-specific codon usage.

#### 2.1 DNA → Protein
- **Function:** Translate a DNA coding sequence into a protein sequence (one-letter or three-letter codes).
- **Options:**
  - **Include Start Codon?** Ensures the first codon is translated as methionine (M).
  - **Include Stop Codon as “*”?** Appends “*” if a stop codon is present.
- **How to use:**
  1. Paste a clean DNA sequence (multiple of 3 bp).
  2. Select output format (one-letter or three-letter amino acids).
  3. Click **Translate**.
  4. View the protein sequence and see which codon table was used.

#### 2.2 Protein → DNA (Reverse Translation)
- **Function:** Reverse-translate a protein sequence into DNA, choosing codons optimized for a selected organism.
- **Options:**
  - **Target Organism:** _E. coli BL21_, _S. cerevisiae_, _P. pastoris_, _H. sapiens_, etc.
  - **GC Preference:** (Optional) Aim for higher or lower GC content.
  - **Avoid Restriction Sites:** Basic list (EcoRI, BamHI, HindIII, XhoI, NdeI).
- **How to use:**
  1. Select **Protein Input** and paste amino acid sequence.
  2. Choose a target organism and any optional constraints.
  3. Click **Reverse Translate**.
  4. Receive a DNA sequence reflecting that organism’s codon usage frequencies.

#### 2.3 Amino Acid Conversions
- **Function:** Convert between one-letter and three-letter amino acid codes.
- **How to use:**
  1. Paste your sequence in one format.
  2. Select target format (One → Three or Three → One).
  3. Click **Convert**.
  4. View the converted sequence with proper delimiters.
"""
        )

    # 3. Codon Optimization
    with st.expander("3. Codon Optimization", expanded=False):
        st.write(
            """
**Purpose:** Optimize DNA sequences for maximal expression in specific host organisms using an advanced algorithm.

#### 3.1 Advanced Algorithm
- **Supported Hosts:**
  - _E. coli BL21_
  - _S. cerevisiae_
  - _P. pastoris_
  - _H. sapiens_
  - _CHO cells_
- **Key Parameters:**
  - **GC Content Target Range:** Slide to choose a minimum and maximum percentage (e.g., 40–60%).
  - **Avoid Restriction Sites:** Multi-select list (EcoRI, BamHI, HindIII, XhoI, NdeI, SalI, XbaI, NotI).
  - **Avoid Repeats (>6 bp):** Check to forbid runs of repeated nucleotides.
  - **Harmonize Usage:** Favor original codons if they already match usage patterns.
- **How to use:**
  1. Paste either a DNA or protein sequence. If you choose **Protein Input**, the tool will translate it first.
  2. Choose your target host, adjust GC range, pick any sites to avoid, and toggle harmonization.
  3. Click **Optimize Sequence**.
  4. **Output Panel** displays:
     - **Optimized Sequence** (full DNA).
     - **Comparison View:** side-by-side codon alignment; changed codons are highlighted.
     - **Stats:** GC% before/after optimization, number of codons changed, verification status (does it translate back correctly?).
     - **Warnings:** Any restriction sites that were successfully avoided, or errors if optimization failed.
"""
        )

    # 4. Extended Synthesis (Long Sequence Fragmentation)
    with st.expander("4. Extended Synthesis (Long Sequence Fragmentation)", expanded=False):
        st.write(
            """
**Purpose:** Fragment long DNA sequences into overlapping oligonucleotides suitable for synthesis and assembly.

#### 4.1 Input Parameters
- **DNA Sequence Input:** Paste directly or upload a FASTA/GenBank/plain-text file.
- **Max Fragment Size (bp):** Choose a size (typical 100–300 bp; default 200).
- **Overlap Size (bp):** Length of overlap between adjacent fragments (default 15).
- **Enzyme Pair:** Select from preset pairs (NdeI/XhoI, BamHI/EcoRI, SalI/XbaI).
- **Cleavage Site (Optional):** Insert protease sites (Thrombin, TEV, Factor Xa, etc.).
- **Advanced Options (under “Advanced Options”):
  - **Show Sequence Warnings:** Flags possible issues (poly-GC runs, long repeats).
  - **Validate Assembly:** Check for GC variation, length variation, and short fragments.

#### 4.2 Running Fragmentation
1. Click **Clean Sequence** to strip invalid characters.
2. Click **Validate** to confirm sequence integrity.
3. Click **Fragment Sequence**. The tool will:
   - Clean and validate your input.
   - Warn if length ≤ fragment size (no fragmentation needed).
   - Ensure overlap < fragment size.
   - Call `fragment_extended_sequence` to generate a list of fragments (each labeled “First,” “Internal,” or “Last”), with:
     - Core sequence (no sticky ends).
     - Forward and reverse primers including sticky-end overhangs and optional cleavage sites.
   - Reassemble overlaps internally to verify they match your original input.

#### 4.3 Viewing Results
- **Metrics Panel** at the top shows total fragments, original length, reassembled length, average fragment size, average GC%.
- **Tabs:**
  - **Fragment Table:** A styled DataFrame listing:
    - Fragment number, type (First/Internal/Last), core length (bp), forward total (bp), reverse total (bp), GC%, snippet of core sequence.  
    - Rows are color-coded: First = green, Internal = blue, Last = red.
  - **Visualization:** Four plots:
    1. Histogram of fragment lengths.
    2. Histogram of fragment GC%.
    3. Bar chart of fragment length vs. fragment number (colors by type).
    4. Line plot of GC% vs. fragment number, with a 50% GC reference line.
  - **Sequences:** Expanders for each fragment show:
    - Core sequence (5′→3′).
    - Forward primer (5′→3′ with sticky end).
    - Reverse primer (5′→3′ with sticky end).
    - GC% and estimated Tₘ (approx. `Tm = 64.9 + 41 × (GC%/100 – 16.4 / length)`).
  - **Export:** Three download buttons:
    - **Download CSV** containing the fragment table.
    - **Download FASTA** with the original sequence plus all forward/reverse primers.
    - **Download Report** (text) summarizing parameters, sequence stats, and fragment details.

"""
        )

    # 5. Hybridization Simulation
    with st.expander("5. Hybridization Simulation", expanded=False):
        st.write(
            """
**Purpose:** Simulate hybridization between two oligos or find internal self-complementarity.

#### 5.1 Alignment
- **Forward Sequence:** Paste your first oligo (or full DNA).
- **Reverse Sequence:** Automatically computed as reverse complement or paste a second oligo.
- **Options:**
  - **Mismatch Penalty:** (Default: –1) Cost for mismatched bases.
  - **Gap Penalty:** (Default: –2) Cost for introducing gaps.
  - **Scoring Matrix:** Simple (+1 match, –1 mismatch) or more advanced.
- **How to use:**
  1. Paste both sequences.
  2. Click **Compute Alignment**.
  3. View:
     - The aligned sequences with “|” marking identical positions.
     - Alignment score and percent identity.
     - Mismatches highlighted in red (if enabled).

#### 5.2 Dot-Matrix Visualization
- **Function:** Creates a dot plot comparing forward vs. reverse sequences.
- **How to use:**
  1. After alignment, click **Show Dot Plot**.
  2. The chart shows diagonal lines where bases match; gaps indicate mismatches.
  3. Hover over the plot (if interactive) to see coordinates.

"""
        )

    # 6. Ligation Check
    with st.expander("6. Ligation Check", expanded=False):
        st.write(
            """
**Purpose:** Quickly verify if two DNA fragments with given enzyme-generated overhangs can ligate.

#### 6.1 Inputs
- **Fragment A Overhang:** Paste or select from enzyme drop-down (e.g., EcoRI “AATT”).
- **Fragment B Overhang:** Paste or select from enzyme drop-down.
- **Molar Ratios/Concentrations:** (Optional) Estimate ligation efficiency.
- **Temperature/Buffer:** (Optional) e.g., T4 ligase conditions at 16°C.

#### 6.2 How to use:
1. Select or paste each 4-nt overhang.
2. Click **Check Compatibility**.
3. View:
   - **Compatibility:** “Compatible” if overhangs are reverse complements; otherwise “Incompatible.”
   - **Predicted ΔG:** Approx. –1.5 kcal/mol per complementary base (–6 kcal for 4 bp), with notes on efficiency.
   - **Schematic:** A simple graphic showing two DNA ends with matching sticky ends.
"""
        )

    # 7. Primer Generator
    with st.expander("7. Primer Generator", expanded=False):
        st.write(
            """
**Purpose:** Design PCR, sequencing, or mutagenesis primers with flexible adapters.

#### 7.1 Inputs
- **Target Sequence:** Paste the template DNA (e.g., plasmid region).
- **Primer Type:** Forward or Reverse.
- **Primer Length:** Desired length (e.g., 18–24 nt).
- **Tm Range:** Minimum and maximum Tₘ (e.g., 58–62°C).
- **GC Content Range:** (Optional) e.g., 40–60%.
- **Prefix/Adapter:** 
  - **Custom 5′ Adapter:** (e.g., “NdeI” adds “CATATG” to the 5′ end).
  - **Special Cases:** If you choose NdeI, the module appends “CATATG” and avoids the “ATG” start codon in the annealing region.

#### 7.2 How to use:
1. Paste or upload a FASTA file containing your template sequence.
2. Choose **Forward** or **Reverse** primer design.
3. Set primer length, Tₘ range, and any GC constraints.
4. (Optional) Select a 5′ adapter from the dropdown (EcoRI, BamHI, HindIII, XhoI, NdeI, etc.).
   - The tool then:
     - Appends the restriction-site sequence to the 5′ end.
     - Recalculates Tₘ including that 5′ extension.
5. Click **Design Primer**.
6. View:
   - **Primer Sequence:** Full 5′ adapter + annealing region.
   - **Annealing Tₘ:** Nearest-neighbor Tₘ for the core region.
   - **GC%:** Overall primer and core region.
   - **Warnings:** Alerts for >4-base G/C runs or predicted hairpins (based on simple self-complementarity).

"""
        )

    # 8. Reverse Complement
    with st.expander("8. Reverse Complement", expanded=False):
        st.write(
            """
**Purpose:** Super-fast utility to clean and reverse complement DNA, with GC%.

#### 8.1 How to use:
1. Paste a DNA sequence (whitespace/newlines are allowed).
2. Click **Clean & Reverse Complement**.
3. View:
   - **Cleaned Sequence:** Uppercase, no spaces/newlines, only A/T/C/G/N.
   - **Reverse Complement:** 5′→3′ string.
   - **GC%:** Percent of G + C in the cleaned sequence.
"""
        )

    # 9. Settings & Session
    with st.expander("9. Settings & Session", expanded=False):
        st.write(
            """
**Purpose:** Customize appearance, manage recent files, and view logs.

#### 9.1 Dark Mode
- **Function:** Switch between light and dark themes.
- **Effect:** Inverts the main background and text color (dark background, light text).

#### 9.2 Recent Files
- **Function:** Tracks recently uploaded FASTA/GenBank files in session state.
- **How it works:**
  - Whenever you upload a file (e.g., in Extended Synthesis), it’s added to **Recent Files**.
  - Click a recent filename to auto-populate the input box.

#### 9.3 Logging
- **Log File:** `g-synth.log` in the working directory.
- **Contents:** 
  - Records each module launch (timestamp + module name).
  - Logs any warnings/errors encountered (e.g., invalid sequence, optimization failure).
"""
        )

    # 10. Ligation Calculator
    with st.expander("10. Ligation Calculator", expanded=False):
        st.write(
            """
**Purpose:** Detailed calculation and visualization of ligation energetics.

#### 10.1 Inputs
- **Insert Overhang:** Drop-down list of enzyme overhangs (EcoRI “AATT,” BamHI “GATC,” etc.) or custom 5′ overhang.
- **Vector Overhang:** Same options as above.
- **Insert:Vector Molar Ratio:** e.g., 3:1, 1:1.
- **Enzyme & Buffer (Optional):** Choose “T4 Ligase – NEB Buffer,” “T7 Ligase – Manufacturer X,” etc.

#### 10.2 How to use:
1. Select or paste each 4-nt overhang.
2. Enter the insert:vector molar ratio (e.g., 3:1).
3. (Optional) Choose enzyme and buffer for a rough efficiency estimate.
4. Click **Calculate Ligase Efficiency**.
5. View:
   - **Compatibility Verdict:** Compatible or Incompatible.
   - **Predicted ΔG:** Using -1.5 kcal/mol per complementary base (e.g., “ΔG_est = -6.0 kcal/mol” for 4 bp).
   - **Relative Efficiency:** A bar chart showing relative efficiency based on ratio and ΔG.
   - **Schematic Diagram:** Simple graphic showing the two overhangs aligned.

"""
        )

    # 11. In Silico Docking
    with st.expander("11. In Silico Docking", expanded=False):
        st.write(
            """
**Purpose:** Basic docking workflow for protein–protein or protein–DNA interactions.

#### 11.1 Structure Prediction
- **Protein Input:**
  - Paste a FASTA or upload a PDB file.
  - If FASTA is provided, uses **ESMFold** (AlphaFold-powered) to generate a model.
- **DNA Input:**
  - Paste a DNA sequence to generate a B-form helix (using 3DNA).

#### 11.2 Docking Engines
- **AutoDock Vina** (Protein–Protein / Protein–Small Molecule):
  - **Binding Site Grid:** Specify coordinates or let the tool guess by sequence motifs.
  - **Parameters:** Exhaustiveness (default 8), number of modes (default 9).
- **HADDOCK-Lite** (Protein–DNA; placeholder):
  - **Restraints:** Option to select binding residues manually.
  - **Outputs:** Docked PDB, interaction score, and cluster ranking.
- **Outputs & Visualization:**
  - **3D Viewer:** Interactive NGL widget (rotatable).
  - **Binding Score:** Vina score (kcal/mol).
  - **Interface Residues:** Table listing contact residues on each partner.
"""
        )

    # 12. Functional Prediction
    with st.expander("12. Functional Prediction", expanded=False):
        st.write(
            """
**Purpose:** Infer protein function using ProtT5 embeddings plus UniProtKB searches and simple classifiers.

#### 12.1 ProtT5 Embedding
- **Input:** Paste a protein FASTA (single chain).
- **Process:** Submits sequence to a ProtT5 model (local or via Hugging Face).
- **Output:** A fixed-length embedding (1024 dimensions) used for downstream tasks.

#### 12.2 UniProtKB Search
- **Function:** Queries the UniProt REST API for similar sequences.
- **Usage:**
  1. After embedding finishes, click **Search UniProtKB**.
  2. View a table of top hits: UniProt ID, percent identity, E-value, and summary of GO annotations.

#### 12.3 Local Classifier (ProtT5-based placeholder)
- **Function:** Uses a pre-trained ProtT5-GO classifier to predict GO Slim categories.
- **Usage:**
  1. Click **Predict Function** after embedding.
  2. See a bar chart of top GO terms with confidence scores.
  3. Click any bar to view details (GO ID, description).
"""
        )

    # 13. Alignment Tools
    with st.expander("13. Alignment Tools", expanded=False):
        st.write(
            """
**Purpose:** Perform pairwise and multiple sequence alignments plus dot-plot visualization.

#### 13.1 Pairwise Alignment
- **Algorithms:**
  - Needleman-Wunsch (global)
  - Smith-Waterman (local)
- **Inputs:**
  - Two sequences (DNA or protein)
  - Scoring matrix:
    - DNA: +1 match, –1 mismatch, –2 gap
    - Protein: BLOSUM62 or PAM250
- **How to use:**
  1. Paste or upload two sequences.
  2. Select algorithm and scoring matrix.
  3. Click **Align**.
  4. View aligned sequences with “|” for identical positions, plus score and percent identity.  
  5. Optionally export alignment as FASTA or CLUSTAL format.

#### 13.2 Multiple Sequence Alignment (MSA)
- **Engines:**
  - MUSCLE (command-line)
  - Clustal Omega (command-line)
- **Inputs:**
  - Upload a FASTA containing at least two sequences.
- **How to use:**
  1. Choose an MSA engine.
  2. (Optional) Set max iterations or output format.
  3. Click **Align Sequences**.
  4. View aligned FASTA and an interactive, color-coded HTML table.  
  5. Download the alignment in FASTA or CLUSTAL format.

#### 13.3 Dot-Plot
- **Function:** Visualize self-similarity or pairwise similarity via a dot matrix.
- **How to use:**
  1. Paste two sequences.
  2. Set window size (e.g., 5 nt) and mismatch threshold.
  3. Click **Generate Dot Plot**.
  4. View a static or interactive dot-plot chart.  
  5. Adjust parameters and regenerate in real time.
"""
        )

    # 14. CRISPR sgRNA Designer
    with st.expander("14. CRISPR sgRNA Designer", expanded=False):
        st.write(
            """
**Purpose:** Identify and score guide RNAs for CRISPR/Cas systems.

#### 14.1 PAM Search
- **Supported Nucleases:**
  - SpCas9 (PAM: NGG)
  - SaCas9 (PAM: NNGRRT)
  - Cas12a/Cpf1 (PAM: TTTV)
  - Cas13 (PAM: NNNN for RNA-targeting)
- **How to use:**
  1. Paste a genomic or gene sequence (DNA).
  2. Select a CRISPR enzyme.
  3. Click **Scan PAM**.
  4. View a table of candidate guides: position, target sequence, PAM sequence, and strand.

#### 14.2 Scoring
- **On-Target (Doench 2016 Score):**  
  - Calculates a score (0–100) for SpCas9 guides.  
  - For other nucleases, uses a placeholder score of 50.
- **Off-Target (Basic):**
  - Runs a 15-nt seed BLAST or Bowtie search locally (placeholder).
  - Reports number of off-target hits with one or fewer mismatches.
- **How to use:**
  1. After PAM scan, select guides.
  2. Click **Score Guides**.
  3. View a table showing each guide’s sequence, PAM, on-target score, off-target hit count, and predicted GC%.

#### 14.3 UCSC Genome Links
- **Function:** Creates a “View in UCSC” link for each guide to open the UCSC Genome Browser at that locus (organism must be specified).
- **How to use:**
  1. Pick a reference genome (e.g., hg38, mm10).
  2. Click **Generate UCSC Links**.
  3. Each guide row includes a clickable link that opens UCSC at the guide’s coordinates.
"""
        )

    # 15. Plasmid Visualizer & Editor
    with st.expander("15. Plasmid Visualizer & Editor", expanded=False):
        st.write(
            """
**Purpose:** Annotate, visualize, and export plasmid maps.

#### 15.1 Input Options
- **GenBank File (.gb/.gbk):**  
  - Upload a GenBank file; the module auto-parses features (genes, CDS, promoters).
- **FASTA + Feature Table:**  
  - Upload a FASTA for your sequence and a CSV specifying feature start, end, type, and label.
- **Raw Sequence:** Paste sequence manually and add features via form.

#### 15.2 Adding, Editing, or Deleting Features
- **Add Feature:**
  1. Click **Add Feature**.
  2. Enter:  
     - Feature type (CDS, Promoter, Terminator, etc.)  
     - Start and end positions  
     - Label  
     - Strand (+/–)  
     - Color (hex code)  
  3. Click **Save**.
- **Edit/Delete Feature:**  
  - Click any feature name in the feature list to edit or delete it.

#### 15.3 Visualization
- **Automatic Mode:**  
  - If sequence length ≥ 2000 bp, show a circular map.  
  - If < 2000 bp, show a linear map.
- **Map Elements:**  
  - Backbone drawn with tick marks every 100 bp.  
  - Features drawn as arrows or boxes, colored by user-chosen colors.  
  - Gene labels placed around the map.
- **Interactive Tools:**  
  - **Zoom In/Out** using on-screen buttons.  
  - **Mouseover Features** to show start/end, length, and annotation.

#### 15.4 Export Options
- **Download:**  
  - **PNG** (raster image)  
  - **SVG** (vector image)  
  - **GenBank** (updated .gbk with any new/edited features)
- **How to use:**
  1. Click your desired download button.  
  2. The file is saved to your local machine.

"""
        )

    st.write("---")
    st.markdown(
        """
Developed by **Dr. Mohamed Merzoug**  
Genomics Technology Platform, Higher School of Biological Sciences of Oran  
Email: **mohamed.merzoug.essbo@gmail.com**  

© 2025 G-Synth Toolkit. All rights reserved.
""",
        unsafe_allow_html=True,
    )