# app/Home.py
#streamlit run app/Home.py --server.address=0.0.0.0 --server.port=8501

import streamlit as st
import os
from utils import save_config
import py3Dmol
from stmol import showmol
import subprocess
from utils import reset_project
import time

st.set_page_config(page_title="OneDock Virtual Screening Pipeline", layout="wide")

# Ensure directories exist
os.makedirs("data/inputs", exist_ok=True)


st.title("Input & Visualization")
st.sidebar.success("Current Step: Input")

# --- A. RECEPTOR SELECTION ---
st.subheader("A. Receptor Input")
structure_known = st.radio("Is the PDB structure of the receptor known?", ["Yes", "No"])

receptor_path = None

if structure_known == "Yes":
    save_config({"structure_known": True})
    c1, c2 = st.columns(2)
    with c1:
        target_file = st.file_uploader("Upload Target Receptor (.pdb)", type="pdb")
        #here we visualize the pdb file
        if target_file is not None:
            pdb_str = target_file.getvalue().decode("utf-8") #decodes bytes into string
            view = py3Dmol.view(width=400, height= 200) #sets the room and its size
            view.addModel(pdb_str, "pdb") #adds protein structure
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo() #zooms to protein
            showmol(view, height = 200, width = 400)
    with c2:
        ref_file = st.file_uploader("Upload Reference Receptor (Optional)", type="pdb")
        #here we visualize the pdb file
        if ref_file is not None:
            pdb_str = ref_file.getvalue().decode("utf-8") #decodes bytes into string
            view = py3Dmol.view(width=400, height= 200) #sets the room and its size
            view.addModel(pdb_str, "pdb") #adds protein structure
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo() #zooms to protein
            showmol(view, height = 200, width = 400)
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

else:
    st.info("Structure Unknown. Generate structures with Bioemu")
    save_config({"structure_known": False})
    uploaded_fasta = st.file_uploader("Upload Target FASTA", type=["fasta"])
    if uploaded_fasta:
        with open("data/inputs/target.fasta", "wb") as f:
            f.write(uploaded_fasta.getbuffer())
        if st.button(" Run BioEmu Prediction"):
            target_output_pdb = "data/inputs/bioemu_target.pdb"
            with st.spinner("BioEmu is generating structures... this may take a minute..."):
                cmd = [
                    "snakemake", 
                    target_output_pdb,
                    "--cores", "1", 
                    "--configfile", "config/config.yaml",
                    "--rerun-incomplete"
                ]

                process = subprocess.run(cmd, capture_output=True, text=True)
                if process.returncode == 0:
                    st.success("Structure generation complete!")
                    save_config({"receptor_path": target_output_pdb})
                else:
                    st.error("BioEmu generation failed.")
                    with st.expander("See Error Log"):
                        st.code(process.stderr)


# --- B. LIGAND LIBRARY ---
os.makedirs("data/inputs/library_split", exist_ok=True)

st.markdown("---")
st.subheader("B. Ligand Library")
smiles_file = st.file_uploader("Upload Ligand Library (.smi / .csv)", type=["smi", "csv"])

if smiles_file:
   # 1. Read the file
    content = smiles_file.getvalue().decode("utf-8")
    lines = [line.strip() for line in content.split('\n') if line.strip()]
    
    st.info(f"Found {len(lines)} ligands. Splitting files...")
    
    # 2. Split into individual files
    for i, line in enumerate(lines):
        smi = line.split()[0]
        fname = f"lig_{i:05d}.smi"
        fpath = os.path.join("data/inputs/library_split", fname)
        
        with open(fpath, "w") as f:
            f.write(smi)
            
    st.success(f"✅ Successfully created {len(lines)} input files in 'data/inputs/library_split/'")

# --- NAVIGATION ---
st.markdown("---")
if st.button("Go to Structure Preparation"):
    st.switch_page("pages/1_Structure_Preparation.py")

st.sidebar.markdown("---")

# Initialize session state for the confirmation button
if "confirm_reset" not in st.session_state:
    st.session_state.confirm_reset = False

def on_reset_click():
    st.session_state.confirm_reset = True

def on_confirm_click():
    status = reset_project()
    if status == "Success":
        st.toast("Project reset successfully!", icon="🗑️")
        time.sleep(1)
        st.rerun() # Refresh the app to clear inputs from memory
    else:
        st.error(status)
    st.session_state.confirm_reset = False

# The Logic
if not st.session_state.confirm_reset:
    st.sidebar.button("🗑️ Reset All Data", on_click=on_reset_click)
else:
    st.sidebar.warning("Are you sure? This will delete all inputs and results.")
    col1, col2 = st.sidebar.columns(2)
    with col1:
        st.button("Yes, Delete", on_click=on_confirm_click, type="primary")
    with col2:
        if st.button("Cancel"):
            st.session_state.confirm_reset = False
            st.rerun()