# app/Home.py
#streamlit run app/Home.py --server.address=0.0.0.0 --server.port=8501

import streamlit as st
import os
import subprocess
import time
import signal
import re
import glob

# --- RDKit Imports ---
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# --- Visualization ---
import py3Dmol
from stmol import showmol

# --- Local Imports ---
from utils import save_config, reset_project, get_available_run_names, load_config

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

            **BioEmu** BioEmu is a generative structural modeling approach used to sample biomolecular conformations. 
            In this project, BioEmu was used only to generate the initial protein structures. No dynamics or free energy 
            calculations were performed with BioEmu; it served exclusively as a tool for initial structure sampling. <br>
            Reference: Lewis, S., Hempel, T., Jiménez-Luna, J., Gastegger, M., Xie, Y., Foong, A. Y., ... & Noé, F. (2025). 
            Scalable emulation of protein equilibrium ensembles with generative deep learning. Science, 389(6761), eadv9817.
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

# --- CONSTANTS (LIGAND PREP) ---
PID_LIG_FILE = "logs/prep/ligprep.pid"
LOG_LIG_FILE = "logs/prep/ligprep.log"

os.makedirs("data/inputs/library_split", exist_ok=True)
os.makedirs("logs/prep", exist_ok=True)


st.markdown("---")
st.subheader("Ligand Library")

# --- 1. CHECK STATUS (LIGAND PREP) ---
lig_job_running = False
if os.path.exists(PID_LIG_FILE):
    try:
        with open(PID_LIG_FILE, "r") as f:
            pid = int(f.read().strip())
        os.kill(pid, 0) # Check if process exists
        lig_job_running = True
    except:
        lig_job_running = False
        # Cleanup dead PID file
        if os.path.exists(PID_LIG_FILE): os.remove(PID_LIG_FILE)

# --- 2. UPLOAD & SPLIT (Only show if not running) ---
if not lig_job_running:
    smiles_file = st.file_uploader("Upload Ligand Library (.smi / .csv)", type=["smi", "csv"])

    if smiles_file:
        # A. GENERATE SECURE PREFIX
        # Get filename without extension (e.g. "My New Lib.smi" -> "My New Lib")
        raw_name = os.path.splitext(smiles_file.name)[0]
        # Replace spaces and weird chars with underscores -> "My_New_Lib"
        library_name = re.sub(r'[^a-zA-Z0-9]', '_', raw_name)

        save_config({"library_name": library_name})
        
        # B. TARGETED CLEANUP
        # Delete only old files belonging to THIS library to prevent "ghost files"
        # e.g., Delete "My_New_Lib_00001.smi" but keep "Old_Lib_00001.smi"
        old_files = glob.glob(f"data/inputs/library_split/{library_name}_*.smi")
        for f in old_files:
            os.remove(f)
            
        if old_files:
            st.info(f"Removed {len(old_files)} old files for library '{library_name}'. Updating...")

        # C. PROCESS CONTENT
        content = smiles_file.getvalue().decode("utf-8")
        lines = content.splitlines()
        
        valid_count = 0
        skipped_count = 0
        skipped_log = []
        progress_bar = st.progress(0)
        
        for i, line in enumerate(lines):
            if not line.strip(): continue
            
            smi = line.split()[0]
            mol = Chem.MolFromSmiles(smi)
            
            # --- FILTER LOGIC ---
            if mol:
                mw = Descriptors.MolWt(mol)
                rot_bonds = Lipinski.NumRotatableBonds(mol)
                
                if rot_bonds > 32:
                    skipped_count += 1
                    skipped_log.append(f"Row {i}: Skipped (>32 rot bonds) - {smi}")
                    continue
                
                if mw > 1000:
                    skipped_count += 1
                    skipped_log.append(f"Row {i}: Skipped (>1000 Da) - {smi}")
                    continue

                # --- SAVE WITH PREFIX ---
                # Format: {LibraryName}_{Number}.smi
                fname = f"{library_name}_{valid_count:05d}.smi" 
                fpath = os.path.join("data/inputs/library_split", fname)
                
                with open(fpath, "w") as f:
                    f.write(smi)
                
                valid_count += 1
            
            if i % 10 == 0:
                progress_bar.progress(min(1.0, i / len(lines)))

        progress_bar.empty()
        
        st.success(f"Saved {valid_count} ligands as '{library_name}_XXXXX.smi'")
        
        if skipped_count > 0:
            st.warning(f"Excluded {skipped_count} ligands.")
            with st.expander("See excluded"):
                st.text("\n".join(skipped_log))

# --- 3. MONITORING MODE (If Job Running) ---
if lig_job_running:
    st.info("Ligand preparation is running in the background...")
    
    # --- NEW: PROGRESS BAR LOGIC ---
    progress_val = 0
    progress_text = "Starting..."
    
    if os.path.exists(LOG_LIG_FILE):
        with open(LOG_LIG_FILE, "r") as f:
            log_content = f.read()
            
            # Find all percentages like "45%" or "100%"
            # pattern: matches one or more digits (\d+) followed by %
            matches = re.findall(r"([\d\.]+)%", log_content)
            
            if matches:
                try:
                    last_pct = float(matches[-1])
                    progress_val = min(last_pct, 100.0)
                    progress_text = f"Docking Progress: {progress_val:.1f}%"
                except ValueError:
                    pass
                
    # Display the bar
    st.progress(progress_val / 100, text=progress_text)
    # -------------------------------

    # Auto-refresh logic
    time.sleep(2) # Wait a bit
    if st.button("Refresh Ligand Status"):
        st.rerun()

    # Show Logs
    log_display = log_content[-2000:] if 'log_content' in locals() else "No logs yet."
            
    with st.expander("Live Preparation Log", expanded=True):
        st.code(log_display)

    # Stop Button
    if st.button("Stop Preparation"):
        with open(PID_LIG_FILE, "r") as f:
            pid = int(f.read().strip())
        try:
            os.kill(pid, signal.SIGTERM)
            st.warning("Preparation job stopped.")
            if os.path.exists(PID_LIG_FILE): os.remove(PID_LIG_FILE)
            st.rerun()
        except Exception as e:
            st.error(f"Could not stop job: {e}")

# --- 4. LAUNCH MODE (If Job Not Running) ---
else:
    # Only show button if files exist
    if os.listdir("data/inputs/library_split"):
        if st.button("Run Ligand Preparation (Background)"):
            
            # Reset logs
            if os.path.exists(LOG_LIG_FILE): os.remove(LOG_LIG_FILE)
            
            with open(LOG_LIG_FILE, "w") as log:
                # Command to run specific rule
                cmd = [
                    "snakemake", 
                    "--cores", "1", 
                    "prepare_all_ligands",  # Your specific target rule
                    "--keep-going",
                    "--rerun-incomplete"
                ]
                
                process = subprocess.Popen(
                    cmd,
                    stdout=log,
                    stderr=log,
                    start_new_session=True # Detach process
                )
                
                # Save PID
                with open(PID_LIG_FILE, "w") as f:
                    f.write(str(process.pid))
            
            st.success("Ligand preparation started! You can scroll down to monitor.")
            st.rerun()

# --- 5. CHECK FOR COMPLETION ---
# If PID file exists but process is dead, it finished
if os.path.exists(PID_LIG_FILE) and not lig_job_running:
    if os.path.exists(LOG_LIG_FILE):
        with open(LOG_LIG_FILE, "r") as f:
            full_log = f.read()
            if "Finished job" in full_log or "Nothing to be done" in full_log:
                st.balloons()
                st.success("Ligand Preparation Complete!")
            else:
                st.error("Ligand Preparation Failed.")
                st.code(full_log[-1000:])
    
    # Cleanup
    if os.path.exists(PID_LIG_FILE): os.remove(PID_LIG_FILE)


#Set run name
current_library_name = load_config().get("library_name")
if current_library_name is not None:
    available_runs = get_available_run_names(current_library_name)

    if not available_runs:
        st.warning("No results found for this library.")
        selected_run = None
    else:
        # Try to select the last used run from config as default
        default_index = 0
        last_run = load_config().get("run_name", "default_run")
        if last_run in available_runs:
            default_index = available_runs.index(last_run)
            
        selected_run = st.selectbox(
            "Select Run / Experiment:", 
            available_runs, 
            index=default_index,
            help="Choose a previous run based on grid parameters or experiment name."
        )
        save_config({"run_name": selected_run})



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
