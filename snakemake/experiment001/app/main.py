import streamlit as st
import os
import yaml
import subprocess
# import py3Dmol  # (Uncomment if you use them later)
# from stmol import showmol
# from pathlib import Path
# import re

# --- SETUP: Ensure ALL directories exist ---
os.makedirs("data/inputs/receptors", exist_ok=True)
os.makedirs("data/inputs/ligands", exist_ok=True) 
os.makedirs("config", exist_ok=True)
os.makedirs("data/results", exist_ok=True)

st.title("🧬 OneDock Pipeline")
st.write("This is a docking pipeline for virtual screening for SLC6A20")

# --- NAVIGATION LOGIC ---
# Initialize the 'step' if it doesn't exist
if 'step' not in st.session_state:
    st.session_state['step'] = 1

# === STEP 1: UPLOAD ===
if st.session_state['step'] == 1:
    st.subheader("Do you have the Receptor Structure available?")
    
    # Using a checkbox to reveal upload
    structure_avail = st.checkbox(label="Known Structure", value=False)

    if structure_avail:
        st.header("Upload Structures")
        
        # Fixed Indentation here
        receptor_file = st.file_uploader(label='Upload your PDB file here', type='pdb')
        
        if receptor_file:
            # Save file to disk
            save_path = os.path.join("data/inputs/receptors", "receptor.pdb")
            with open(save_path, "wb") as f:
                f.write(receptor_file.getbuffer())
            st.success(f"Receptor saved to {save_path}!")

            # The Navigation Button
            # Only show this if file is uploaded
            if st.button("Go to Structure Preparation ->"):
                st.session_state['step'] = 2
                st.rerun() # Forces the app to reload and show Step 2

# === STEP 2: PREPARATION ===
elif st.session_state['step'] == 2:
    st.header("Structure Preparation")
    st.write("Now configuring the preparation parameters...")
    
    # Add your preparation widgets here
    clean_pdb = st.checkbox("Remove Water Molecules?", value=True)
    
    if st.button("<- Back"):
        st.session_state['step'] = 1
        st.rerun()
        
    if st.button("Run Preparation Pipeline"):
        st.write("Running Snakemake...")
        # subprocess.run(...)
    


# # --- SECTION 1: SIDEBAR (Configuration) ---
# st.sidebar.header("1. Grid Configuration")
# st.sidebar.info("Define the search space for Vina.")

# # Create a dictionary to hold config values
# config_data = {}
# config_data['center_x'] = st.sidebar.number_input("Center X", value=0.0)
# config_data['center_y'] = st.sidebar.number_input("Center Y", value=0.0)
# config_data['center_z'] = st.sidebar.number_input("Center Z", value=0.0)
# st.sidebar.markdown("---")
# config_data['size_x'] = st.sidebar.number_input("Size X", value=20.0)
# config_data['size_y'] = st.sidebar.number_input("Size Y", value=20.0)
# config_data['size_z'] = st.sidebar.number_input("Size Z", value=20.0)
# config_data['exhaustiveness'] = st.sidebar.slider("Exhaustiveness", 8, 32, 8)

# # --- SECTION 2: FILE UPLOADS ---
# st.header("2. Upload Structures")
# c1, c2 = st.columns(2)

# with c1:
#     receptor_file = st.file_uploader(label='Upload your PDB file here', type='pdb')
#     if receptor_file:
#         # Save file to disk so Snakemake can see it
#         with open(os.path.join("data/inputs", "receptor.pdb"), "wb") as f:
#             f.write(receptor_file.getbuffer())
#         st.success("Receptor saved!")

# with c2:
#     ligand_file = st.file_uploader("Upload your Ligand Library here", type="smi")
#     if ligand_file:
#         with open(os.path.join("data/inputs", "ligand.smi"), "wb") as f:
#             f.write(ligand_file.getbuffer())
#         st.success("Ligand saved!")

# pocket_file = st.file_uploader(label='Upload your pocket PDB file here', type='pdb')

# # --- SECTION 3: EXECUTION ---
# st.header("3. Run Docking")

# if st.button("🚀 Start Docking Pipeline"):
#     if not (receptor_file and ligand_file):
#         st.error("Please upload both files first.")
#     else:
#         # A. Save the configuration to YAML
#         with open("config/vina_config.yaml", "w") as f:
#             yaml.dump(config_data, f)
        
#         # B. Run Snakemake
#         with st.spinner("Running AutoDock Vina... (This may take a minute)"):
#             # Note: We use --cores 1 and point to our Snakefile
#             # We assume your Snakefile is at the root or workflow/Snakefile
#             cmd = [
#                 "snakemake",
#                 "--cores", "1",
#                 "--configfile", "config/vina_config.yaml",
#                 # Add --use-conda or --use-singularity here if needed
#             ]
            
#             process = subprocess.run(cmd, capture_output=True, text=True)
            
#             if process.returncode == 0:
#                 st.balloons()
#                 st.success("Docking Finished Successfully!")
                
#                 # Check if results exist and show them
#                 if os.path.exists("data/results/output.pdbqt"):
#                     st.text("Preview of output:")
#                     with open("data/results/output.pdbqt", "r") as res:
#                         st.code(res.read(1000) + "...") # Show first 1000 chars
#             else:
#                 st.error("Pipeline Failed")
#                 st.text("Error Log:")
#                 st.code(process.stderr)