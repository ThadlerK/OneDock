import streamlit as st
import pandas as pd
import os

st.title("Results & Ranking")

result_file = "data/results/docking_report.csv"

if os.path.exists(result_file):
    df = pd.read_csv(result_file)
    st.subheader("Docking Scores")
    st.dataframe(df.sort_values(by="Affinity_kcal_mol"))
    
    st.subheader("Best Pose Visualization")
    st.write("Select a ligand to visualize (Placeholder logic)")
    
    # --- NAVIGATION TO ADME ---
    st.markdown("---")
    st.info("Proceed to ADME screening to evaluate drug-likeness properties of top candidates")
    if st.button("Go to ADME Screening →"):
        st.switch_page("pages/ADME_Screening.py")
else:
    st.warning("No results found. Please run the docking pipeline first.")