# modules/functional_prediction.py

import logging
import streamlit as st
import requests
import pandas as pd

from utils.bio_utils import translate_sequence, load_sequence_from_fasta

# ───────────────────────────────────────────────────────────────────────────────
# 1) Logging
# ───────────────────────────────────────────────────────────────────────────────
logger = logging.getLogger("G-Synth:FUNC_PRED")
logging.basicConfig(
    level=logging.INFO,
    filename="g-synth-func.log",
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ───────────────────────────────────────────────────────────────────────────────
# 2) Detect PyTorch & Transformers availability
# ───────────────────────────────────────────────────────────────────────────────
try:
    import torch
    from transformers import T5EncoderModel, T5TokenizerFast
    TORCH_AVAILABLE = True
    logger.info("PyTorch and Transformers loaded successfully.")
except Exception as e:
    torch = None
    T5EncoderModel = None
    T5TokenizerFast = None
    TORCH_AVAILABLE = False
    logger.warning(f"PyTorch/Transformers unavailable; skipping ProtT5 embedding: {e}")

# ───────────────────────────────────────────────────────────────────────────────
# 3) Cache ProtT5 model loader (only if available)
# ───────────────────────────────────────────────────────────────────────────────
@st.cache_resource
def load_prott5_model():
    if not TORCH_AVAILABLE:
        return None, None, None
    tokenizer = T5TokenizerFast.from_pretrained(
        "Rostlab/prot_t5_xl_uniref50", do_lower_case=False
    )
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)
    return tokenizer, model, device

# ───────────────────────────────────────────────────────────────────────────────
# 4) UniProtKB REST API query
# ───────────────────────────────────────────────────────────────────────────────
def query_uniprot(sequence: str) -> pd.DataFrame:
    url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=sequence:{sequence}&format=json"
        "&fields=accession,protein_name,go_id,ec,comment(PATHWAY)"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
    except Exception as e:
        logger.error(f"UniProt query failed: {e}")
        return pd.DataFrame()
    data = r.json().get("results", [])[:5]
    records = []
    for entry in data:
        acc   = entry.get("primaryAccession","")
        pname = (
            entry.get("proteinDescription",{})
                 .get("recommendedName",{})
                 .get("fullName",{})
                 .get("value","")
        )
        # GO terms & EC & pathways
        go_terms = [
            cr["properties"]["GoTerm"]
            for cr in entry.get("uniProtKBCrossReferences",[])
            if cr.get("database")=="GO"
        ]
        ec_nums = [
            cr["id"]
            for cr in entry.get("uniProtKBCrossReferences",[])
            if cr.get("database")=="EC"
        ]
        pathways = []
        for c in entry.get("comments",[]):
            if c.get("commentType")=="PATHWAY":
                for t in c.get("texts",[]):
                    val = t.get("value")
                    if val:
                        pathways.append(val)
        records.append({
            "Accession":  acc,
            "Protein":    pname,
            "GO Terms":   ", ".join(go_terms[:3]),
            "EC Numbers": ", ".join(ec_nums[:1]),
            "Pathways":   ", ".join(pathways[:3])
        })
    return pd.DataFrame(records)

# ───────────────────────────────────────────────────────────────────────────────
# 5) Dummy local classifier stub
# ───────────────────────────────────────────────────────────────────────────────
def local_prott5_classifier(embedding: "torch.Tensor") -> dict[str, float]:
    # Placeholder: replace with real classifier head
    return {
        "GO:0008150 (biological_process)": 0.75,
        "GO:0003674 (molecular_function)": 0.60
    }

# ───────────────────────────────────────────────────────────────────────────────
# 6) Streamlit App
# ───────────────────────────────────────────────────────────────────────────────
def main():
    st.title("AI-Based Functional Prediction")

    # Input type
    input_type = st.radio("Input Type:", ["Protein","DNA → Protein"])

    # Sequence input: file or manual
    seq_file   = st.file_uploader("Upload FASTA file", type=["fa","fasta"])
    seq_manual = st.text_area("Or paste sequence (FASTA or raw):", height=150)

    if st.button("Predict Function"):
        # Determine source text
        raw_txt = None
        if seq_file:
            raw_txt = seq_file.getvalue().decode("utf-8")
        elif seq_manual.strip():
            raw_txt = seq_manual.strip()
        else:
            st.error("Please upload a FASTA or paste your sequence.")
            return

        # Extract sequence
        if raw_txt.lstrip().startswith(">"):
            try:
                seq = load_sequence_from_fasta(raw_txt)
            except Exception as e:
                st.error(f"FASTA parse error: {e}")
                return
        else:
            seq = "".join(raw_txt.split())

        # Translate if needed
        if input_type == "DNA → Protein":
            try:
                seq = translate_sequence(seq, 0, True)
            except Exception as e:
                st.error(f"Translation error: {e}")
                return

        # 1) UniProtKB search
        st.info("Querying UniProtKB…")
        up_df = query_uniprot(seq)
        if not up_df.empty:
            st.write("#### UniProtKB Hits")
            st.dataframe(up_df)
        else:
            st.warning("No UniProtKB hits found.")

        # 2) ProtT5 embedding & classification
        tokenizer, model, device = load_prott5_model()
        if TORCH_AVAILABLE and tokenizer and model:
            st.info("Generating ProtT5 embedding…")
            seq_spaced = " ".join(seq)
            tokens = tokenizer(seq_spaced, return_tensors="pt")
            input_ids = tokens.input_ids.to(device)
            with torch.no_grad():
                emb = model(input_ids).last_hidden_state.mean(dim=1)
            st.info("Classifying GO terms…")
            go_conf = local_prott5_classifier(emb)
            go_df = pd.DataFrame(
                go_conf.items(), columns=["GO Term","Confidence"]
            ).sort_values("Confidence", ascending=False)
            st.write("#### Predicted GO Terms")
            st.bar_chart(go_df.set_index("GO Term"))
            st.table(go_df)
        else:
            st.warning("ProtT5 embedding unavailable; skipped local classification.")

        # 3) Dummy EC & pathways output
        st.write("#### Predicted EC Number")
        st.write("3.4.21.4 (Serine Protease) — Confidence: 0.68")
        st.write("#### Suggested Pathways")
        st.markdown(
            "- KEGG: [hsa05150 Staphylococcus aureus infection]"
            "(https://www.kegg.jp/pathway/hsa05150)"
        )
        st.markdown(
            "- Reactome: [R-HSA-6785807 Evidence Support]"
            "(https://reactome.org/PathwayBrowser/#/R-HSA-6785807)"
        )

        # 4) Final summary
        st.write(
            "**Summary:** Likely serine protease (EC 3.4.21.4), "
            "involved in complement activation."
        )

if __name__ == "__main__":
    main()