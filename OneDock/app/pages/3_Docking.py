import streamlit as st
import subprocess
import yaml
from utils import save_config
from utils import load_config
import shutil
import os

st.title("Run Docking")

# --- ENGINE SELECTION ---
engine = st.selectbox("Select Docking Engine", ["AutoDock Vina", "Schrödinger Glide"])

if engine == "Schrödinger Glide":
    st.warning("⚠️ Ensure you have a valid Schrödinger license in the Docker container.")

save_config({"docking_engine": engine})

# --- EXECUTION ---
st.markdown("---")
st.subheader("Pipeline Execution")

if st.button("Launch Docking Pipeline"):
    try:
        with open("config/config.yaml", "r") as f:
            config = yaml.safe_load(f)
            
        pocket_residues = config.get("pocket_residues") 

        if not pocket_residues or not str(pocket_residues).strip():
            st.error("**Stop:** No binding pocket defined!")
            st.info("Please go to the **Structure** page and select residues (or enter them manually) before docking.")
            st.stop()  
            
    except FileNotFoundError:
        st.error("Config file not found. Please save settings first.")
        st.stop()

    # 2. Spinner & Cleanup (Passiert nur, wenn oben nicht gestoppt wurde)
    with st.spinner(f"Running docking pipeline..."):
        
        # Cleanup der alten Ergebnisse
        if os.path.exists("data/results/stats"):
            shutil.rmtree("data/results/stats")
        if os.path.exists("data/results/logs"):
            shutil.rmtree("data/results/logs")
        if os.path.exists("data/results/poses"):
            shutil.rmtree("data/results/poses")
        if os.path.exists("data/results/docking_report.csv"):
            os.remove("data/results/docking_report.csv")
        
        # 3. Construct Command
        cmd = [
            "snakemake",
            "data/results/docking_report.csv",
            "--cores", "1",
            "--configfile", "config/config.yaml",
            "--rerun-incomplete"
        ]
        
        # 4. Run
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        # 5. Handle Result
        if process.returncode == 0:
            st.balloons()
            st.success("Docking Completed Successfully!")
            st.switch_page("pages/4_Results.py")
        else:
            st.error("Pipeline Failed.")
            with st.expander("See Error Log"):
                st.code(process.stderr)