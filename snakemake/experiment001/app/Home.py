# app/Home.py
import streamlit as st
import os
from utils import save_config
import py3Dmol
from stmol import showmol
import subprocess

st.set_page_config(page_title="OneDock Virtual Screening Pipeline", layout="wide")

# Ensure directories exist
os.makedirs("data/inputs", exist_ok=True)

st.title("1️. Input & Visualization")
st.sidebar.success("Current Step: Input")

# --- A. RECEPTOR SELECTION ---
st.subheader("A. Receptor Input")
structure_known = st.radio("Is the PDB structure of the receptor known?", ["Yes", "No"])

receptor_path = None

if structure_known == "Yes":
    c1, c2 = st.columns(2)
    with c1:
        target_file = st.file_uploader("Upload Target Receptor (.pdb)", type="pdb")
    with c2:
        ref_file = st.file_uploader("Upload Reference Receptor (Optional)", type="pdb")

    if target_file:
        # Save to disk
        receptor_path = os.path.join("data/inputs", "target.pdb")
        with open(receptor_path, "wb") as f:
            f.write(target_file.getbuffer())
        
        st.success(f"Saved: {receptor_path}")
        save_config({"receptor_path": receptor_path})

        # --- VISUALIZATION (Only if uploaded) ---
        # st.write("### 🧬 Structure Preview")
        # pdb_str = target_file.getvalue().decode("utf-8")
        # view = py3Dmol.view(width=800, height=400)
        # view.addModel(pdb_str, "pdb")
        # view.setStyle({'cartoon': {'color': 'spectrum'}})
        # view.zoomTo()
        # showmol(view, height=400, width=800)

        if st.button("Run Preparation"):
            with st.spinner("Running Python Preparation Script..."):
                target_output_file = "data/interim/target.pdbqt"
                
                cmd = ["snakemake", "--cores", "1", target_output_file]
                process = subprocess.run(cmd, capture_output=True, text=True)
                
                if process.returncode == 0:
                    st.success("Preparation complete!")
                    # Visualization logic here...
                else:
                    st.error("Preparation failed.")
                    st.code(process.stderr)

else:
    st.warning("Structure Unknown")
    if st.button("🧬 Run BioEmu Prediction"):
        st.info("Triggering BioEmu... (Placeholder)")
        # In real life: subprocess.run(["snakemake", "bioemu_target"])

# --- B. LIGAND LIBRARY ---
st.markdown("---")
st.subheader("B. Ligand Library")
smiles_file = st.file_uploader("Upload Ligand Library (.smi / .csv)", type=["smi", "csv"])

if smiles_file:
    lib_path = os.path.join("data/inputs", "library.smi")
    with open(lib_path, "wb") as f:
        f.write(smiles_file.getbuffer())
    save_config({"ligand_library": lib_path})
    st.success("Library Saved.")

# --- NAVIGATION ---
st.markdown("---")
if st.button("Go to Structure Preparation ➡️"):
    st.switch_page("pages/1_Structure_Preparation.py")