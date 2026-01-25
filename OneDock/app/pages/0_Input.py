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

st.set_page_config(page_title="Input", layout="wide")

# Ensure directories exist
os.makedirs("data/inputs", exist_ok=True)


st.title("Input & Visualization")

# --- A. RECEPTOR SELECTION ---
st.subheader("Receptor Input")
st.write("Upload your receptor and ligand library. If you want you can \n" \
        "add a reference receptor. Make sure your receptors are clean and \n" \
        "don't contain any ligands, water molecules or additional proteins.")

structure_known = st.radio("Is the PDB structure of the receptor known?", ["Yes", "No"])

receptor_path = None

if structure_known == "Yes":
    save_config({"structure_known": True})
    c1, c2 = st.columns(2)
    with c1:
        target_file = st.file_uploader("Upload Target Receptor (.pdb)", type="pdb")
        # Save to disk
        target_path = os.path.join("data/inputs", "target.pdb")
        if target_file:
            with open(target_path, "wb") as f:
                f.write(target_file.getbuffer())
            save_config({"target_path": target_path})

        #here we visualize the pdb file
        if os.path.isfile(target_path):
            with open(target_path, "r") as f:
                target_str = f.read()

            view = py3Dmol.view(width=400, height= 200) #sets the room and its size
            view.addModel(target_str, "pdb") #adds protein structure
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo() #zooms to protein
            showmol(view, height = 200, width = 400)

            if st.button("Run Receptor Preparation"):
                with st.spinner("Running Python Preparation Script..."):
                    target_output_file = "data/interim/target_prep.pdbqt"
                    
                    cmd = ["snakemake", "--cores", "1", target_output_file, "--rerun-incomplete"]
                    process = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if process.returncode == 0:
                        st.success("Preparation complete!")
                        # Visualization logic here...
                    else:
                        st.error("Preparation failed.")
                        st.code(process.stderr)

    with c2:
        ref_file = st.file_uploader("Upload Reference Receptor (Optional)", type="pdb")
        # Save to disk
        ref_path = os.path.join("data/inputs", "reference.pdb")
        if ref_file:
            with open(ref_path, "wb") as f:
                f.write(ref_file.getbuffer())
            save_config({"ref_path": ref_path})

        #here we visualize the pdb file
        if os.path.isfile(ref_path):
            with open(ref_path, "r") as f:
                ref_str = f.read()

            view = py3Dmol.view(width=400, height= 200) #sets the room and its size
            view.addModel(ref_str, "pdb") #adds protein structure
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo() #zooms to protein
            showmol(view, height = 200, width = 400)

            if st.button("Run Reference Preparation"):
                with st.spinner("Running Python Preparation Script..."):
                    target_output_file = "data/interim/reference_prep.pdbqt"
                    
                    cmd = ["snakemake", "--cores", "1", target_output_file, "--rerun-incomplete"]
                    process = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if process.returncode == 0:
                        st.success("Preparation complete!")
                    else:
                        st.error("Preparation failed.")
                        st.code(process.stderr)


else:
    st.info("""
            Since you don't know the structure of your receptor, we will use Bio Emu to
            predict it.\n

            **BioEmu** is a tool...
            """)
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
st.subheader("Ligand Library")
smiles_file = st.file_uploader("Upload Ligand Library (.smi / .csv)", type=["smi", "csv"])

if smiles_file:
   # 1. Read the file
    content = smiles_file.getvalue().decode("utf-8")
    lines = content.splitlines() 
    lines = [line.strip() for line in lines if line.strip()]
    
    
    # 2. Split into individual files
    for i, line in enumerate(lines):
        smi = line.split()[0]
        fname = f"lig_{i:05d}.smi"
        fpath = os.path.join("data/inputs/library_split", fname)
        
        with open(fpath, "w") as f:
            f.write(smi)


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
    st.switch_page("pages/1_Pocket_Detection.py")



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