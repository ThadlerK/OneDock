import streamlit as st
import subprocess
import yaml
from utils import save_config

st.title("Run Docking")

# --- ENGINE SELECTION ---
engine = st.selectbox("Select Docking Engine", ["AutoDock Vina", "Schrödinger Glide"])

if engine == "Schrödinger Glide":
    st.warning("⚠️ Ensure you have a valid Schrödinger license in the Docker container.")

save_config({"docking_engine": engine})

# --- EXECUTION ---
st.markdown("---")
st.subheader("Pipeline Execution")

if st.button("🚀 Launch Docking Pipeline"):
    with st.spinner(f"Running {engine}..."):
        # 1. Construct Command
        # We point Snakemake to our config and the desired output file
        cmd = [
            "snakemake",
            "data/results/docking_report.csv",
            "--cores", "1",
            "--configfile", "config/config.yaml"
        ]
        
        # 2. Run
        process = subprocess.run(cmd, capture_output=True, text=True)
        
        # 3. Handle Result
        if process.returncode == 0:
            st.balloons()
            st.success("Docking Completed Successfully!")
            st.switch_page("pages/4_Results.py")
        else:
            st.error("Pipeline Failed.")
            with st.expander("See Error Log"):
                st.code(process.stderr)