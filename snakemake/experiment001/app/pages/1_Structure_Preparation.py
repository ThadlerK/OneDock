import streamlit as st
from utils import save_config
import subprocess

st.title("Structure Preparation")
st.subheader("Receptor conversion")

if st.button("Run Receptor Preparation"):
    with st.spinner("Running Python Preparation Script..."):
        target_output_file = "data/interim/receptor_ready.pdbqt"
        
        cmd = ["snakemake", "--cores", "1", target_output_file]
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            st.success("Preparation complete!")
            # Visualization logic here...
        else:
            st.error("Preparation failed.")
            st.code(process.stderr)

st.subheader("Ligand conversion")

if st.button("Run Ligand Preparation"):
    with st.spinner("Converting all ligands to PDBQT..."):

        target_rule = "prepare_all_ligands"
        
        # We pass the rule name directly to snakemake
        cmd = ["snakemake", "--cores", "1", target_rule]
        
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        if process.returncode == 0:
            st.success("All ligands converted successfully!")
        else:
            st.error("Ligand Preparation failed.")
            with st.expander("Error Log"):
                st.code(process.stderr)

# --- NAVIGATION ---
st.markdown("---")
if st.button("Go to Pocket Definition step"):
    st.switch_page("pages/2_Pocket_Prediction.py")