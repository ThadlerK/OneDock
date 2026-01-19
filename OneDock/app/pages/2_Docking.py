import streamlit as st
import subprocess
import yaml
import shutil
import os
import signal
import time

# --- CONSTANTS ---
PID_FILE = "docking.pid"
LOG_FILE = "docking.log"
STATUS_FILE = "docking_status.txt" # "running", "success", "failed"

st.title("Run Docking")

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
    # Read the final log to guess status (or use Snakemake exit codes if we wrapped it)
    log_tail = get_log_content()
    cleanup_job() # Reset for next run
    
    if "Finished job 0." in log_tail or "100%" in log_tail or "Nothing to be done" in log_tail:
        st.balloons()
        st.success(" Docking Job Completed while you were away!")
        st.switch_page("pages/3_Results.py")
    else:
        st.error("The background job failed.")
        with st.expander("Error Log"):
            st.code(log_tail)

# 2. DISPLAY ACTIVE JOB (MONITORING MODE)
if job_running:
    st.info("Docking pipeline is running in the background...")
    st.caption("You can close this tab or VS Code. The job will continue.")
    
    # Auto-refresh the page every 10 seconds to check status
    st.empty() # Placeholder
    time.sleep(1) # Small delay
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
    
    if st.button("Launch Docking Pipeline"):
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
        if os.path.exists("data/results"):
            shutil.rmtree("data/results")
            os.makedirs("data/results", exist_ok=True)
            
        # --- C. Launch Background Process ---
        # We redirect stdout/stderr to a file so we can read it later
        with open(LOG_FILE, "w") as log:
            # Construct the command
            cmd = [
                "snakemake",
                "--cores", "1",
                "--configfile", "config/config.yaml",
                "--rerun-incomplete",
                "--keep-going" # Don't stop on first error
            ]
            
            # Popen starts the process without waiting for it
            process = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=log,
                start_new_session=True # CRITICAL: Detaches from the parent (Streamlit)
            )
            
            # Save the PID
            with open(PID_FILE, "w") as f:
                f.write(str(process.pid))
                
        st.success("Job started in background!")
        st.rerun() # Reload page to switch to "Monitoring Mode"