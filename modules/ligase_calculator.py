# modules/ligase_calculator.py

"""
Ligation Calculator for G-Synth

Step 1: Calculates required insert DNA mass at common molar ratios.
Step 2: Calculates reaction setup volumes given final volume and DNA concentrations.
"""

import streamlit as st
import pandas as pd

# ───────────────────────────────────────────────────────────────────────────────
# UNIT CONVERSION FACTORS
# ───────────────────────────────────────────────────────────────────────────────
LENGTH_UNITS = {
    "bp": 1,
    "kb": 1e3,
    "Mb": 1e6
}
MASS_UNITS = {
    "ng": 1,
    "µg": 1e3
}

# ───────────────────────────────────────────────────────────────────────────────
# RENDER FUNCTION
# ───────────────────────────────────────────────────────────────────────────────
def render_ligation_calculator():
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
        <h1 style="margin: 0; font-size: 2.5rem; font-weight: 800;">🔗 Ligation Calculator</h1>
        <p style="margin: 0.5rem 0 0 0; font-size: 1.2rem; opacity: 0.9;">
            Professional DNA ligation planning with molar ratio calculations
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown(
        "Calculate the required insert DNA mass at several molar ratios, then "
        "compute the reaction setup volumes given the final reaction volume and "
        "DNA concentrations."
    )

    # ── STEP 1: Required insert mass ──────────────────────────────────────────────
    st.subheader("### 📊 Step 1: Insert Mass Calculation")
    col1, col2, col3 = st.columns(3)
    with col1:
        ins_len  = st.number_input(
            "Insert DNA length", min_value=0.0, value=2.0, step=0.1, format="%.2f"
        )
        ins_unit = st.selectbox("", list(LENGTH_UNITS.keys()), index=1, key="ins_len_unit")
    with col2:
        vec_len  = st.number_input(
            "Vector DNA length", min_value=0.0, value=6.0, step=0.1, format="%.2f"
        )
        vec_unit = st.selectbox("", list(LENGTH_UNITS.keys()), index=1, key="vec_len_unit")
    with col3:
        vec_mass = st.number_input(
            "Vector DNA mass", min_value=0.0, value=100.0, step=0.1, format="%.2f"
        )
        mass_unit = st.selectbox("", list(MASS_UNITS.keys()), index=0, key="vec_mass_unit")

    # Convert to base units
    ins_bp      = ins_len  * LENGTH_UNITS[ins_unit]
    vec_bp      = vec_len  * LENGTH_UNITS[vec_unit]
    vec_mass_ng = vec_mass * MASS_UNITS[mass_unit]

    # Compute insert mass at standard ratios
    ratios = [1, 2, 3, 5, 7]
    rows = []
    for r in ratios:
        if vec_bp > 0:
            ins_mass = vec_mass_ng * (r * ins_bp / vec_bp)
            ins_mass_str = f"{ins_mass:.2f}"
        else:
            ins_mass_str = "—"
        rows.append({
            "Ratio": f"{r}:1",
            "Insert mass (ng)": ins_mass_str
        })

    df_mass = pd.DataFrame(rows).set_index("Ratio")
    st.table(df_mass)

    st.markdown(
        "<details>\n"
        "<summary>🔎 Formula</summary>\n\n"
        "**Insert mass (ng)** = (insert:vector ratio) × (vector mass (ng)) × "
        "(insert length / vector length)\n"
        "</details>",
        unsafe_allow_html=True
    )

    # ── STEP 2: Reaction Volume Calculation ──────────────────────────────────────
    st.subheader("### ⚗️ Step 2: Reaction Setup & Volumes")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        choice = st.selectbox(
            "Choose ratio", df_mass.index.tolist(), index=2, key="ratio_select"
        )
        r_val = float(choice.split(":")[0])
    with c2:
        final_vol = st.number_input(
            "Final reaction volume (µL)", min_value=1.0, value=20.0, step=0.5
        )
    with c3:
        vec_conc = st.number_input(
            "Vector conc. (ng/µL)", min_value=0.1, value=50.0, step=0.1
        )
    with c4:
        ins_conc = st.number_input(
            "Insert conc. (ng/µL)", min_value=0.1, value=50.0, step=0.1
        )

    enzyme_vol = st.number_input(
        "T4 DNA Ligase (µL)", min_value=0.0, value=1.0, step=0.1
    )

    # Calculate volumes
    ins_mass = vec_mass_ng * (r_val * ins_bp / vec_bp) if vec_bp > 0 else 0.0
    vec_vol  = vec_mass_ng / vec_conc  if vec_conc > 0 else 0.0
    ins_vol  = ins_mass      / ins_conc if ins_conc > 0 else 0.0
    buffer_vol = final_vol / 10.0
    water_vol  = final_vol - (vec_vol + ins_vol + buffer_vol + enzyme_vol)

    if water_vol < 0:
        st.error("🚨 Volumes exceed total reaction volume—please adjust inputs.")
    else:
        df_setup = pd.DataFrame({
            "Component":   ["Vector", "Insert", "10× Buffer", "T4 DNA Ligase", "Water"],
            "Volume (µL)": [
                round(vec_vol, 2),
                round(ins_vol, 2),
                round(buffer_vol, 2),
                round(enzyme_vol, 2),
                round(water_vol, 2)
            ]
        })
        st.table(df_setup)
        
        # Enhanced download with blue styling
        csv = df_setup.to_csv(index=False)
        st.download_button(
            label="📥 Download reaction setup (CSV)",
            data=csv,
            file_name="ligation_setup.csv",
            mime="text/csv",
            type="primary"
        )

# ───────────────────────────────────────────────────────────────────────────────
# ENTRYPOINT FUNCTIONS (REQUIRED)
# ───────────────────────────────────────────────────────────────────────────────

def main():
    """Main entrypoint function for the Ligation Calculator module."""
    render_ligation_calculator()

def app():
    """Alternative entrypoint function for compatibility."""
    render_ligation_calculator()

# For standalone execution
if __name__ == "__main__":
    main()
