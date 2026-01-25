import streamlit as st
import subprocess
import yaml
import shutil
import os
import signal
import time
import re
from utils import save_config

st.set_page_config(page_title = "Docking")

# --- CONSTANTS ---
PID_FILE = "logs/docking/docking.pid"
LOG_FILE = "logs/docking/docking.log"
STATUS_FILE = "docking_status.txt" # "running", "success", "failed"

st.title("Docking with AutoDock Vina")
st.info("""
        **AutoDock Vina** is a fast molecular docking tool that predicts binding poses and 
        binding affinities of small molecules to protein targets. It explores a user-defined 
        binding site, generates multiple ligand poses, and scores them using an empirical 
        scoring function. \n

        Reference: Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy 
        of docking with a new scoring function, efficient optimization, and multithreading. 
        Journal of computational chemistry, 31(2), 455-46
        """)

# --- HELPER FUNCTIONS ---

def is_job_running():
    """Checks if the PID in the pid file is still active."""
    if not os.path.exists(PID_FILE):
        return False
    
    try:
        with open(PID_FILE, "r") as f:
            pid = int(f.read().strip())
        
        # Signal 0 does not kill the process, just checks if it exists
        os.kill(pid, 0)
        return True
    except (OSError, ValueError):
        # Process is dead or file is corrupt
        return False

def get_log_content():
    """Reads the last few lines of the log file."""
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, "r") as f:
            return f.read()[-2000:] # Last 2000 chars
    return "No logs yet."

def cleanup_job():
    """Removes tracking files."""
    if os.path.exists(PID_FILE):
        os.remove(PID_FILE)
    if os.path.exists(STATUS_FILE):
        os.remove(STATUS_FILE)

# --- UI LOGIC ---

# 1. CHECK CURRENT STATUS
# If the app reloads, we check if a job was already running
job_running = is_job_running()

# If the job was running but the PID is gone now, it finished (success or fail)
if os.path.exists(PID_FILE) and not job_running:
    # Read the final log to guess status
    log_tail = get_log_content()
    cleanup_job() # Reset for next run
    
    if "Finished job 0." in log_tail or "100%" in log_tail or "Nothing to be done" in log_tail:
        st.balloons()
        st.success(" Docking Job Completed while you were away!")
        st.switch_page("pages/3_Docking Results.py")
    else:
        st.error("The background job failed.")
        with st.expander("Error Log"):
            st.code(log_tail)

# 2. DISPLAY ACTIVE JOB (MONITORING MODE)
if job_running:
    st.info("Docking pipeline is running in the background...")
    st.caption("You can close this tab or VS Code. The job will continue.")
    
    # --- NEW: PROGRESS BAR LOGIC ---
    progress_val = 0
    progress_text = "Initializing docking..."
    
    # We read the log to find the percentage
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, "r") as f:
            log_content = f.read()
            # Regex: find all digits followed by '%' (e.g., 50%)
            matches = re.findall(r"([\d\.]+)%", log_content)
            
            if matches:
                # Take the last occurrence and convert to float
                try:
                    last_pct = float(matches[-1])
                    progress_val = min(last_pct, 100.0)
                    progress_text = f"Docking Progress: {progress_val:.1f}%"
                except ValueError:
                    pass
    
    # Render the bar
    st.progress(progress_val / 100, text=progress_text)
    # -------------------------------
    
    # Auto-refresh the page every 2 seconds
    st.empty() 
    time.sleep(2) 
    if st.button("Refresh Status"):
        st.rerun()

    with st.expander("Live Log Output", expanded=True):
        st.code(get_log_content())
        
    if st.button("Stop Job"):
        with open(PID_FILE, "r") as f:
            pid = int(f.read().strip())
        try:
            os.kill(pid, signal.SIGTERM)
            st.warning("Job stopped.")
            cleanup_job()
            st.rerun()
        except Exception as e:
            st.error(f"Could not stop job: {e}")

# 3. DISPLAY LAUNCH BUTTON (INPUT MODE)
else:
    st.subheader("Pipeline Execution")

    input_dir = "data/inputs/library_split"
    available_libs = set()
    
    if os.path.exists(input_dir):
        files = os.listdir(input_dir)
        for f in files:
            if f.endswith(".smi"):
                # Extract prefix (e.g. "Chembl" from "Chembl_0001.smi")
                # We assume format is Name_Number.smi
                match = re.match(r"(.+)_\d+\.smi", f)
                if match:
                    available_libs.add(match.group(1))
    
    # Sort and create dropdown
    lib_list = sorted(list(available_libs))
    
    if not lib_list:
        st.error("No ligand libraries found! Please go to Ligand Prep page.")
        st.stop()
        
    # Get current selection from config if it exists
    current_lib = lib_list[0]
    if os.path.exists("config/config.yaml"):
        with open("config/config.yaml") as f:
            c = yaml.safe_load(f) or {}
            if c.get("library_name") in lib_list:
                current_lib = c.get("library_name")

    selected_lib = st.selectbox(
        "Select Ligand Library to Dock", 
        lib_list, 
        index=lib_list.index(current_lib)
    )
    
    st.info(f"Selected library: **{selected_lib}**")
    
    # --- A. NEW: DOCKING SETTINGS ---
    
    if os.path.exists("config/config.yaml"):
        with open("config/config.yaml", "r") as f:
            existing_conf = yaml.safe_load(f) or {}
            current_grid = existing_conf.get("grid_size", 20)
            current_exhaust = existing_conf.get("exhaustiveness", 8)

    with st.expander("Docking Parameters", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            # Box Size Input
            grid_size = st.number_input(
                "Box Size (Å)", 
                min_value=10, 
                max_value=30, 
                value=int(current_grid),
                help="Size of the cubic box centered on the pocket (default: 20Å)"
            )
        with col2:
            # Exhaustiveness Input
            exhaustiveness = st.number_input(
                "Exhaustiveness", 
                min_value=1, 
                max_value=100, 
                value=int(current_exhaust),
                help="How hard Vina searches. Higher = slower but more accurate (default: 8)"
            )
            
    # --------------------------------
    
    if st.button("Launch Docking Pipeline"):
        save_config({
            "grid_size": grid_size,
            "exhaustiveness": exhaustiveness
        })

        # --- A. Validation ---
        try:
            with open("config/config.yaml", "r") as f:
                config = yaml.safe_load(f)
            
            pocket_residues = config.get("pocket_residues") 
            if not pocket_residues or not str(pocket_residues).strip():
                st.error("**Stop:** No binding pocket defined!")
                st.stop()
                
            ref_path = config.get("ref_path")
            ref_residues = config.get("ref_residues")
            if ref_path and (not ref_residues or not str(ref_residues).strip()):
                 st.error("**Stop:** Reference path set but no reference residues defined!")
                 st.stop()
                 
        except FileNotFoundError:
            st.error("Config file not found.")
            st.stop()

        # --- B. Cleanup Old Results ---
        # Note: Be careful here. If you want to keep old results (different libraries),
        # you might want to skip this delete or make it more targeted.
        # if os.path.exists("data/results"):
        #     shutil.rmtree("data/results")
        #     os.makedirs("data/results", exist_ok=True)
            
        # --- C. Launch Background Process ---
        with open(LOG_FILE, "w") as log:
            cmd = [
                "snakemake",
                "--cores", "1",
                "--configfile", "config/config.yaml",
                "--rerun-incomplete",
                "--keep-going"
            ]
            
            process = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=log,
                start_new_session=True 
            )
            
            with open(PID_FILE, "w") as f:
                f.write(str(process.pid))
                
        st.success("Job started in background!")
        st.rerun()