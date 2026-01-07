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
else:
    st.warning("No results found. Please run the docking pipeline first.")