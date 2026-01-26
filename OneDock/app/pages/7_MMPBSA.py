import streamlit as st
import subprocess
import yaml
import os
import glob
import pandas as pd
import signal
import time

st.set_page_config(page_title="MMPBSA Analysis", layout="wide")
st.title("MMPBSA Rescoring")

os.makedirs("logs", exist_ok=True)

# --- CONSTANTS ---
PID_FILE = "logs/MMPBSA/mmpbsa.pid"
LOG_FILE = "logs/MMPBSA/mmpbsa.log"
CONFIG_FILE = "config/config.yaml"

# --- HELPER FUNCTIONS ---
def is_job_running():
    """Checks if the background process is still alive."""
    if not os.path.exists(PID_FILE):
        return False
    try:
        with open(PID_FILE, "r") as f:
            pid = int(f.read().strip())
        os.kill(pid, 0) # Signal 0 checks existence
        return True
    except (OSError, ValueError):
        return False

def get_log_content():
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, "r") as f:
            return f.read()[-5000:] # Show last 5000 chars
    return "No logs yet."

def stop_job():
    if os.path.exists(PID_FILE):
        try:
            with open(PID_FILE, "r") as f:
                pid = int(f.read().strip())
            os.kill(pid, signal.SIGTERM)
        except:
            pass
        os.remove(PID_FILE)

# --- 1. LOAD DATA ---
try:
    with open(CONFIG_FILE) as f:
        config = yaml.safe_load(f)
except:
    config = {}

lib_name = config.get("library_name", "target")
RUN_NAME = config.get("run_name", 'default_run')
DOCKING_RESULTS = f"data/results/docking_report_target_{lib_name}_{RUN_NAME}.csv"

df = None # Initialisieren mit None
if os.path.exists(DOCKING_RESULTS):
    df = pd.read_csv(DOCKING_RESULTS)
else:
    st.warning(f"No docking results found for library '{lib_name}' ({DOCKING_RESULTS}).")

# --- 2. CHECK JOB STATUS ---
job_active = is_job_running()

# Detect if job just finished (File exists but process dead)
if os.path.exists(PID_FILE) and not job_active:
    os.remove(PID_FILE) # Cleanup
    st.success("Background job finished!")
    # Optional: You could parse logs here to see if it actually succeeded

# --- 3. UI LOGIC ---

if job_active:
    # --- MONITORING MODE ---
    st.info("MMPBSA Pipeline is running in the background...")
    st.caption("You can close this tab. The calculations will continue on the server.")
    
    # Auto-refresh button
    if st.button("Refresh Logs"):
        st.rerun()
        
    with st.expander("Live Log Output", expanded=True):
        st.code(get_log_content())
        
    if st.button("Stop Job"):
        stop_job()
        st.warning("Job stopped.")
        st.rerun()

else:
    # --- INPUT MODE ---
    if df is not None:
        
        # 1. DELTA RANK BERECHNUNG (falls Referenz vorhanden)
        # Wir versuchen, die Datei basierend auf der Config zu finden
        ref_lib_name = config.get("library_name", "")
        REF_FILE = f"data/results/docking_report_reference_{lib_name}_{RUN_NAME}.csv"
        
        has_rank_metric = False
        
        if os.path.exists(REF_FILE):
            try:
                df_r = pd.read_csv(REF_FILE)
                # Sicherstellen, dass Spaltennamen passen
                if "Affinity_kcal_mol" in df_r.columns:
                    df_r = df_r.rename(columns={"Affinity_kcal_mol": "Affinity_Ref"})
                
                # Ränge berechnen (Niedrigere Energie = Besserer Rang -> 1 ist bester)
                df["Rank_Target"] = df["Affinity_kcal_mol"].rank(method="min", ascending=True)
                df_r["Rank_Ref"] = df_r["Affinity_Ref"].rank(method="min", ascending=True)
                
                # Mergen
                df = pd.merge(df, df_r[["Ligand", "Affinity_Ref", "Rank_Ref"]], on="Ligand", how="left")
                
                # Delta Rank berechnen (Positiv = Verbesserung im Target)
                # Bsp: Ref Rank 100 (schlecht) - Target Rank 10 (gut) = +90 Gain
                df["Delta_Rank"] = df["Rank_Ref"] - df["Rank_Target"]
                has_rank_metric = True
            except Exception as e:
                st.warning(f"Konnte Referenzdaten nicht verarbeiten: {e}")

        # 2. UI EINGABEN
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Selection Strategy")
            
            # Auswahl: Top N
            default_top_n = config.get("mmpbsa", {}).get("top_n", 5)
            top_n = st.number_input("Auto-Select Top N Ligands", min_value=0, max_value=50, value=default_top_n, help="Sorting by **Highest Rank Gain** (Delta Rank)")
            
            # Info, wonach sortiert wird
            if has_rank_metric:
                st.success(f"Sorting by **Highest Rank Gain** (Delta Rank)")
            else:
                st.info("Sorting by **Affinity** (Reference not found)")

            # ZUSÄTZLICH: Manuelle Auswahl
            all_ligands = sorted(df["Ligand"].unique().tolist())
            manual_selection = st.multiselect(
                "Add specific ligands manually (optional)", 
                options=all_ligands,
                help="These ligands are addtionally analysed."
            )

        with col2:
            st.subheader("MD Parameters")
            ranks_per_ligand = st.slider("Poses per Ligand", 1, 5, 1)
            heat_steps = st.number_input("Heat Steps", value=config.get("mmpbsa", {}).get("md_steps_heat", 20000))
            prod_steps = st.number_input("Production Steps", value=config.get("mmpbsa", {}).get("md_steps_prod", 100000))

        # 3. LISTE ERSTELLEN
        # A. Automatische Auswahl (Top N)
        if top_n > 0:
            if has_rank_metric:
                # Sortieren nach Delta Rank (Absteigend -> Höchster Gewinn zuerst)
                # NaNs (keine Ref-Daten) werden ans Ende geschoben
                df_sorted = df.sort_values("Delta_Rank", ascending=False)
            else:
                # Fallback: Sortieren nach Affinity (Aufsteigend -> Niedrigste Energie zuerst)
                df_sorted = df.sort_values("Affinity_kcal_mol", ascending=True)
            
            auto_ligands = df_sorted.head(top_n)["Ligand"].tolist()
        else:
            auto_ligands = []

        # B. Zusammenführen (Auto + Manuell) und Duplikate entfernen
        final_selection = list(set(auto_ligands + manual_selection))
        
        # Anzeigen
        if final_selection:
            st.write(f"**Selected Ligands ({len(final_selection)}):** {', '.join(final_selection)}")
            
            # Warnung bei großen Jobs
            total_runs = len(final_selection) * ranks_per_ligand
            if total_runs > 20:
                st.warning(f"⚠️ High workload: {total_runs} simulations. This might take a long time.")
            else:
                st.info(f"Total simulations: {len(final_selection)} ligands × {ranks_per_ligand} poses = {total_runs} runs")
        else:
            st.warning("No ligands selected via Top N or Manual input.")

        if st.button("Launch MMPBSA Pipeline (Background)"):
            # 1. Save Config
            if "mmpbsa" not in config: config["mmpbsa"] = {}
            config["mmpbsa"]["top_n"] = top_n
            config["mmpbsa"]["md_steps_heat"] = heat_steps
            config["mmpbsa"]["md_steps_prod"] = prod_steps
            with open(CONFIG_FILE, "w") as f:
                yaml.dump(config, f)

            # 2. Build Targets List
            targets = []
            for lig in top_ligands:
                for r in range(1, ranks_per_ligand + 1):
                    targets.append(f"data/results/mmpbsa/{lig}_rank{r}/FINAL_RESULTS_MMPBSA.dat")

            # 3. Launch Process
            with open(LOG_FILE, "w") as log:
                cmd = ["snakemake", "--cores", "1", "--use-conda", "--keep-going", "--rerun-incomplete"] + targets
                
                process = subprocess.Popen(
                    cmd,
                    stdout=log,
                    stderr=log,
                    start_new_session=True # Detach process
                )
                
                # Save PID
                with open(PID_FILE, "w") as f:
                    f.write(str(process.pid))
                    
            st.success("Job started! Page will reload...")
            time.sleep(1)
            st.rerun()

# --- 4. RESULTS PREVIEW ---
st.divider()
st.subheader("Results Table")
results_files = glob.glob("data/results/mmpbsa/*/FINAL_RESULTS_MMPBSA.dat")

if results_files:
    data = []
    for f in results_files:
        try:
            with open(f) as txt:
                content = txt.read()
                for line in content.splitlines():
                    if line.startswith("DELTA TOTAL"):
                        val = float(line.split()[2])
                        parts = f.split("/")
                        # Extract folder name: {lig}_rank{r}
                        folder = parts[-2]
                        lig_id = folder.split("_rank")[0]
                        rank_id = folder.split("_rank")[1]
                        
                        data.append({
                            "Ligand": lig_id, 
                            "Rank": rank_id, 
                            "Delta G (MMPBSA)": val
                        })
                        break
        except:
            pass
            
    if data:
        df_res = pd.DataFrame(data).sort_values("Delta G (MMPBSA)")
        st.dataframe(df_res, use_container_width=True)
        
        # Simple Download
        csv = df_res.to_csv(index=False).encode('utf-8')
        st.download_button("Download MMPBSA results as CSV", csv, "mmpbsa_results.csv", "text/csv")
    else:
        st.info("Parsing results found, but could not extract Delta G values yet.")
else:
    st.caption("No MMPBSA results generated yet.")