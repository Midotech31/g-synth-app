# modules/plasmid_visualizer.py

"""
Plasmid Visualizer & Editor
- Inputs: GenBank file upload or raw sequence + feature-table CSV/TSV
- Manual “Add Feature” form (name, start, end, strand, color)
- Parse GenBank via Biopython, maintain list of SeqFeature objects
- On “Export”: write new GenBank
- Visualization: dna_features_viewer.CircularGraphicRecord or GraphicRecord
- Download PNG/SVG
"""

import streamlit as st
import tempfile
import os
from io import BytesIO
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import CircularGraphicRecord, GraphicRecord, GraphicFeature
import pandas as pd
import logging

logger = logging.getLogger("G-Synth:PLASMID_VIS")

def parse_genbank(file_content: str):
    """
    Parse GenBank content and return SeqRecord and list of features.
    """
    record = SeqIO.read(file_content, "genbank")
    return record, record.features

def main():
    st.title("Plasmid Visualizer & Editor")

    gb_upload = st.file_uploader("Upload GenBank File:", type=["gb","gbk"])
    seq_input = st.text_area("Or Paste Raw Sequence (FASTA):", height=100)
    feat_input = st.file_uploader("Upload Feature-Table CSV/TSV (optional):", type=["csv","tsv"])

    record = None
    features = []

    if gb_upload:
        temp = tempfile.NamedTemporaryFile(delete=False)
        temp.write(gb_upload.getbuffer())
        temp.flush()
        record, features = parse_genbank(temp.name)
    elif seq_input:
        seq = "".join(seq_input.split())
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        record = SeqRecord(Seq(seq), id="plasmid1", name="plasmid1", description="")
        if feat_input:
            df = pd.read_csv(feat_input) if feat_input.name.endswith(".csv") else pd.read_table(feat_input)
            for _, row in df.iterrows():
                start = int(row["Start"])
                end = int(row["End"])
                name = row["Feature Name"]
                strand = 1 if row["Strand"] == "+" else -1
                color = row["Color"]
                feat = SeqFeature(FeatureLocation(start, end), type="misc_feature", qualifiers={"label": name, "color": color})
                features.append(feat)
    else:
        st.info("Upload a GenBank file or paste raw sequence.")
        return

    # Display existing features
    st.write("### Current Features")
    feat_list = []
    for idx, f in enumerate(features, 1):
        label = f.qualifiers.get("label", [f.type])[0]
        start = int(f.location.start)
        end = int(f.location.end)
        strand = "+" if f.location.strand == 1 else "-"
        color = f.qualifiers.get("color", ["#777777"])[0]
        feat_list.append({"No.": idx, "Name": label, "Start": start, "End": end, "Strand": strand, "Color": color})
    if feat_list:
        st.table(pd.DataFrame(feat_list))

    # Add new feature form
    st.write("---")
    st.subheader("Add / Edit Feature")
    name = st.text_input("Feature Name:")
    start = st.number_input("Start (0-based):", min_value=0, value=0, step=1)
    end = st.number_input("End (1-based):", min_value=1, value=1, step=1)
    strand = st.radio("Strand:", ["+", "-"])
    color = st.color_picker("Color (hex):", "#FF0000")
    if st.button("Add Feature"):
        feat = SeqFeature(FeatureLocation(start, end), type="misc_feature", qualifiers={"label": name, "color": color})
        features.append(feat)
        st.success("Feature added. Rerun to see updated table.")

    # Visualization
    st.write("---")
    st.subheader("Plasmid Map")
    length = len(record.seq)
    if length < 2000:
        # Linear visualization
        graphic_feats = []
        for f in features:
            lab = f.qualifiers.get("label", [f.type])[0]
            col = f.qualifiers.get("color", ["#777777"])[0]
            graphic_feats.append(GraphicFeature(start=int(f.location.start), end=int(f.location.end), strand=int(f.location.strand), color=col, label=lab))
        graphic_record = GraphicRecord(sequence_length=length, features=graphic_feats)
        ax, _ = graphic_record.plot(figure_width=10)
        buf = BytesIO()
        ax.figure.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        st.image(buf)
        # Download PNG
        st.download_button("Download PNG", buf, file_name="plasmid_linear.png", mime="image/png")
    else:
        # Circular visualization
        graphic_feats = []
        for f in features:
            lab = f.qualifiers.get("label", [f.type])[0]
            col = f.qualifiers.get("color", ["#777777"])[0]
            graphic_feats.append(GraphicFeature(start=int(f.location.start), end=int(f.location.end), strand=int(f.location.strand), color=col, label=lab))
        circ_record = CircularGraphicRecord(sequence_length=length, features=graphic_feats)
        ax, _ = circ_record.plot(figure_width=8)
        buf = BytesIO()
        ax.figure.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        st.image(buf)
        # Download PNG
        st.download_button("Download PNG", buf, file_name="plasmid_circular.png", mime="image/png")

    # Export new GenBank
    if st.button("Export GenBank"):
        out_path = os.path.join(tempfile.gettempdir(), "edited_plasmid.gb")
        record.features = features
        SeqIO.write(record, out_path, "genbank")
        with open(out_path, "rb") as f:
            st.download_button("Download GenBank", f, file_name="edited_plasmid.gb")
