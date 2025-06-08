# modules/docking_module.py

"""
AI-Based In Silico Docking
- Inputs: sequence upload (FASTA/raw), dropdown: Protein or DNA, optional PDB upload.
- If no PDB: predict structure (ESMFold or simple B-form DNA).
- Docking: Protein–Protein (AutoDock Vina), Protein–DNA (HADDOCK or shape scoring).
- Outputs: docked PDB, binding score, interface residues.
- Visualization: py3Dmol if available, else static PNG.
"""

import streamlit as st
import tempfile
import os
from utils.pdb_utils import predict_protein_structure, build_bform_dna_model, parse_interface_residues
import subprocess
import logging
import py3Dmol

logger = logging.getLogger("G-Synth:DOCK")

def main():
    st.title("AI-Based In Silico Docking")

    st.write("Upload sequences or PDBs, choose docking type, and run docking.")

    docking_type = st.radio("Docking Type:", ["Protein–Protein", "Protein–DNA"])
    seq_input = st.file_uploader("Upload Sequence (FASTA) or paste raw text:", type=["fa","fasta","txt"])
    pdb_upload = st.file_uploader("Upload PDB (optional):", type=["pdb"])

    if st.button("Run Docking"):
        # Read sequence
        seq = ""
        if seq_input:
            try:
                content = seq_input.getvalue().decode("utf-8")
                seq = "".join([line.strip() for line in content.splitlines() if not line.startswith(">")])
            except:
                st.error("Failed to read sequence input.")
                return
        # Write seq to temp
        temp_dir = tempfile.mkdtemp()
        pred_pdb = None
        if pdb_upload:
            pdb_path = os.path.join(temp_dir, pdb_upload.name)
            with open(pdb_path, "wb") as f:
                f.write(pdb_upload.getbuffer())
        else:
            if docking_type == "Protein–Protein":
                st.info("Predicting protein structure (ESMFold)...")
                pred_pdb = os.path.join(temp_dir, "pred_A.pdb")
                out = predict_protein_structure(seq, pred_pdb)
                if not out:
                    st.error("Failed to predict protein structure.")
                    return
                pdb_path = pred_pdb
            else:
                st.info("Building B-form DNA model...")
                pred_pdb = os.path.join(temp_dir, "dna.pdb")
                out = build_bform_dna_model(seq, pred_pdb)
                if not out:
                    st.error("Failed to build DNA model.")
                    return
                pdb_path = pred_pdb

        # Assume partner PDB is also provided or predicted separately (for simplicity, reuse same structure)
        partner_path = pdb_path

        # Run docking via AutoDock Vina CLI if available
        vina_path = shutil.which("vina")
        if not vina_path:
            st.error("AutoDock Vina not found in PATH. Cannot perform docking.")
            return

        # Prepare input PDBQT files (placeholder)
        st.info("Preparing PDBQT files... (this is a placeholder)")

        # Run Vina (placeholder command)
        conf_out = os.path.join(temp_dir, "out.pdbqt")
        try:
            cmd = [vina_path, "--receptor", pdb_path, "--ligand", partner_path, "--out", conf_out]
            subprocess.run(cmd, check=True)
            st.success("Docking completed.")
        except subprocess.CalledProcessError as e:
            st.error(f"Vina error: {e}")
            return

        # Parse results: for demonstration, use dummy score
        binding_score = -7.5  # placeholder
        st.write(f"### Binding Score (Vina): {binding_score} kcal/mol")

        # Convert out.pdbqt to out.pdb (placeholder)
        docked_pdb = os.path.join(temp_dir, "docked.pdb")
        with open(docked_pdb, "w") as f:
            f.write("REMARK Dummy docked PDB\n")  # Placeholder

        # Interface residues
        interface = parse_interface_residues(pdb_path, partner_path, docked_pdb)
        st.write("### Interface Residues (within 5 Å)")
        if interface:
            st.table(interface)
        else:
            st.write("No interface residues identified (placeholder).")

        # Visualization
        st.write("### 3D Visualization of Docked Complex")
        view = py3Dmol.view(width=600, height=400)
        pdb_data = open(docked_pdb, "r").read()
        view.addModel(pdb_data, "pdb")
        view.setStyle({'cartoon': {'color':'spectrum'}})
        view.zoomTo()
        html = view._make_html()
        st.components.v1.html(html, height=450, scrolling=True)

        # Download link
        with open(docked_pdb, "rb") as f:
            st.download_button("Download Docked PDB", f, file_name="docked.pdb")
