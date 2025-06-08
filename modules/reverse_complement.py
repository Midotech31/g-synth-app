# reverse_complement.py

import streamlit as st
import re

from utils.bio_utils import calculate_gc, calculate_tm_consensus

def clean_sequence(seq: str) -> str:
    """
    Remove any characters not A/T/C/G/N and convert to uppercase.
    """
    return re.sub(r"[^ATCGN]", "", seq.strip().upper())

def complement(seq: str) -> str:
    """
    Return the complement of a DNA sequence (A<->T, C<->G, N<->N).
    Assumes the input is already cleaned (uppercase, only A/T/C/G/N).
    """
    trans_table = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(trans_table)

def reverse(seq: str) -> str:
    """
    Return the reverse of the DNA sequence.
    """
    return seq[::-1]

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    return complement(reverse(seq))

def main():
    """
    Entry point for the Reverse Complement module.
    Should be called by app.py when this module is selected.
    """

    # Title and description (no st.set_page_config here)
    st.title("DNA Sequence Tools: Reverse / Complement / Reverse Complement")
    st.write(
        """
        Use the tabs below to compute the **Reverse**, **Complement**, or **Reverse Complement**
        of any DNA sequence. Invalid characters (anything not A/T/C/G/N) will be removed automatically.
        """
    )

    tabs = st.tabs(["Reverse", "Complement", "Reverse Complement"])

    # ─────────────────────────────────────────────────────────────────────────
    # Tab 1: Reverse
    # ─────────────────────────────────────────────────────────────────────────
    with tabs[0]:
        st.header("Reverse Sequence")
        st.write("Enter a DNA sequence and click **Reverse** to obtain its reversed string.")

        with st.form(key="form_reverse", clear_on_submit=True):
            seq_input = st.text_area(
                label="Input DNA Sequence (5'-3')",
                placeholder="e.g., ATGCGTACGTTAGC",
                height=120,
                key="input_rev",
            )
            submit_rev = st.form_submit_button("Reverse")

        if submit_rev:
            cleaned = clean_sequence(seq_input)
            if not cleaned:
                st.error("Please enter a valid DNA sequence (A/T/C/G/N).")
            else:
                rev_seq = reverse(cleaned)
                st.subheader("Reversed Sequence (5'-3'):")
                st.code(rev_seq, language="plaintext")
                st.markdown(f"**Length:** {len(rev_seq)} bp")


    # ─────────────────────────────────────────────────────────────────────────
    # Tab 2: Complement
    # ─────────────────────────────────────────────────────────────────────────
    with tabs[1]:
        st.header("Complement Sequence")
        st.write("Enter a DNA sequence and click **Complement** to obtain its complement (3'-5').")

        with st.form(key="form_complement", clear_on_submit=True):
            seq_input = st.text_area(
                label="Input DNA Sequence (5'-3')",
                placeholder="e.g., ATGCGTACGTTAGC",
                height=120,
                key="input_comp",
            )
            submit_comp = st.form_submit_button("Complement")

        if submit_comp:
            cleaned = clean_sequence(seq_input)
            if not cleaned:
                st.error("Please enter a valid DNA sequence (A/T/C/G/N).")
            else:
                comp_seq = complement(cleaned)
                gc_content = calculate_gc(comp_seq)
                st.subheader("Complement Sequence (3'-5'):")
                st.code(comp_seq, language="plaintext")
                st.markdown(
                    f"**Length:** {len(comp_seq)} bp  \n"
                    f"**GC Content:** {gc_content:.1f}%"
                )


    # ─────────────────────────────────────────────────────────────────────────
    # Tab 3: Reverse Complement
    # ─────────────────────────────────────────────────────────────────────────
    with tabs[2]:
        st.header("Reverse Complement Sequence")
        st.write(
            "Enter a DNA sequence and click **Reverse Complement** to obtain its reverse complement, "
            "along with GC content and estimated melting temperature (Tm)."
        )

        with st.form(key="form_revcomp", clear_on_submit=True):
            seq_input = st.text_area(
                label="Input DNA Sequence (5'-3')",
                placeholder="e.g., ATGCGTACGTTAGC",
                height=120,
                key="input_revcomp",
            )
            submit_revcomp = st.form_submit_button("Reverse Complement")

        if submit_revcomp:
            cleaned = clean_sequence(seq_input)
            if not cleaned:
                st.error("Please enter a valid DNA sequence (A/T/C/G/N).")
            else:
                revcomp_seq = reverse_complement(cleaned)
                gc_content = calculate_gc(revcomp_seq)
                tm_estimate = calculate_tm_consensus(revcomp_seq)

                st.subheader("Reverse Complement (5'-3'):")
                st.code(revcomp_seq, language="plaintext")
                st.markdown(
                    f"**Length:** {len(revcomp_seq)} bp  \n"
                    f"**GC Content:** {gc_content:.1f}%  \n"
                    f"**Estimated Tm:** {tm_estimate:.1f} C"
                )