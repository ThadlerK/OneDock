import streamlit as st
import pandas as pd
import altair as alt
import os
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, JsCode
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(layout="wide", page_title="Docking Results")

# --- KONFIGURATION ---
TARGET_FILE = "data/results/docking_report_target.csv"
REF_FILE = "data/results/docking_report_reference.csv"

st.title("Results")

# --- 1. DATEN LADEN & MERGEN ---
@st.cache_data
def load_data(target_path, ref_path):
    if not os.path.exists(target_path):
        return None
    
    # Target laden
    df_t = pd.read_csv(target_path)
    df_t = df_t.rename(columns={"Affinity_kcal_mol": "Affinity_Target"})
    df_t["Rank_Target"] = df_t["Affinity_Target"].rank(method="min", ascending=True)

    # Referenz laden (falls vorhanden)
    if os.path.exists(ref_path):
        df_r = pd.read_csv(ref_path)
        df_r = df_r.rename(columns={"Affinity_kcal_mol": "Affinity_Ref"})
        df_r["Rank_Ref"] = df_r["Affinity_Ref"].rank(method="min", ascending=True)
        
        # Merge auf Ligand Name
        merged = pd.merge(df_t, df_r[["Ligand", "Affinity_Ref", "Rank_Ref"]], on="Ligand", how="inner")
        
        # Metriken berechnen
        merged["Delta_Affinity"] = merged["Affinity_Target"] - merged["Affinity_Ref"]
        merged["Delta_Rank"] = merged["Rank_Ref"] - merged["Rank_Target"]
        
        # Sortieren nach bester Spezifität (Delta)
        return merged.sort_values(by="Delta_Affinity")
    else:
        # Nur Target -> Sortieren nach Affinität
        return df_t.sort_values(by="Affinity_Target")

df = load_data(TARGET_FILE, REF_FILE)

if df is None:
    st.warning("No docking results found. Please run the pipeline first.")
    st.stop()

# --- 2. ANALYSE & PLOTS ---
has_reference = "Affinity_Ref" in df.columns

# Wir nutzen einen Expander für die gesamte Analyse
with st.expander("Analysis Charts", expanded=True):
    
    # Auswahl-Element für den Plot-Typ
    plot_type = st.radio(
        "Select Visualization:", 
        ["Specificity Scatter Plots", "Affinity Distributions"], 
        horizontal=True
    )
    
    st.divider()

    # --- OPTION A: SCATTER PLOTS (Wie vorher) ---
    if plot_type == "Specificity Scatter Plots":
        if has_reference:
            col_chart1, col_chart2 = st.columns(2)
            
            with col_chart1:
                st.markdown("**Affinity Comparison** (Lower Right = Better Specificity)")
                chart_scatter = alt.Chart(df).mark_circle(size=60).encode(
                    x=alt.X('Affinity_Ref', title='Ref Affinity (kcal/mol)'),
                    y=alt.Y('Affinity_Target', title='Target Affinity (kcal/mol)'),
                    color=alt.Color('Delta_Affinity', scale=alt.Scale(scheme='redblue', reverse=True), title='Specificity'),
                    tooltip=['Ligand', 'Affinity_Target', 'Affinity_Ref', 'Delta_Affinity']
                ).interactive()
                st.altair_chart(chart_scatter, use_container_width=True)

            with col_chart2:
                st.markdown("**Rank Comparison** (Top Left = Rank Improvement)")
                chart_rank = alt.Chart(df).mark_circle(size=60).encode(
                    x=alt.X('Rank_Ref', title='Rank on Reference'),
                    y=alt.Y('Rank_Target', title='Rank on Target'),
                    color=alt.Color('Delta_Rank', scale=alt.Scale(scheme='viridis'), title='Rank Gain'),
                    tooltip=['Ligand', 'Rank_Target', 'Rank_Ref', 'Delta_Rank']
                ).interactive()
                st.altair_chart(chart_rank, use_container_width=True)
        else:
            st.info("Scatter plots require a reference receptor to calculate specificity.")

    # --- OPTION B: HISTOGRAMME (Neu) ---
    elif plot_type == "Affinity Distributions":
        
        # 1. Daten "schmelzen" (Wide to Long Format) für Altair
        # Altair mag es, wenn alle Werte in einer Spalte stehen und eine zweite Spalte die Kategorie angibt
        cols_to_melt = ["Affinity_Target"]
        if has_reference:
            cols_to_melt.append("Affinity_Ref")
            
        df_melted = df.melt(
            id_vars=["Ligand"], 
            value_vars=cols_to_melt, 
            var_name="Receptor", 
            value_name="Energy (kcal/mol)"
        )
        
        # Namen für die Legende aufhübschen
        df_melted["Receptor"] = df_melted["Receptor"].replace({
            "Affinity_Target": "Target Receptor", 
            "Affinity_Ref": "Reference Receptor"
        })

        # 2. Histogramm erstellen
        base = alt.Chart(df_melted)
        
        hist = base.mark_bar(opacity=0.6, binSpacing=0).encode(
            alt.X("Energy (kcal/mol)", bin=alt.Bin(maxbins=40), title="Binding Affinity (kcal/mol)"),
            alt.Y("count()", title="Number of Ligands", stack=None), # stack=None lässt die Balken überlappen
            alt.Color("Receptor", scale=alt.Scale(scheme="category10")),
            tooltip=["Receptor", "count()"]
        ).properties(
            height=350,
            title="Distribution of Binding Affinities"
        ).interactive()

        st.altair_chart(hist, use_container_width=True)
        
        # Statistische Zusammenfassung anzeigen
        st.markdown("#### Summary Statistics")
        stats_cols = st.columns(2 if has_reference else 1)
        
        with stats_cols[0]:
            st.caption("Target Receptor")
            st.write(df["Affinity_Target"].describe().to_frame().T)
            
        if has_reference:
            with stats_cols[1]:
                st.caption("Reference Receptor")
                st.write(df["Affinity_Ref"].describe().to_frame().T)

st.divider()

# --- 3. INTERAKTIVE TABELLE & STRUKTUR ---
st.subheader("Interactive Ranking & Visualization")
st.caption("Click on a row in the table to see the chemical structure.")

# Spalten auswählen, die wir anzeigen wollen
if has_reference:
    display_cols = ["Ligand", "Affinity_Target", "Affinity_Ref", "Delta_Affinity", "Delta_Rank", "Smiles"]
else:
    display_cols = ["Ligand", "Affinity_Target", "Rank_Target", "Smiles"]

df_display = df[display_cols].copy()

# Layout: 2/3 Tabelle, 1/3 Bild
col_table, col_viz = st.columns([2, 1])

with col_table:
    # AgGrid Konfiguration
    gb = GridOptionsBuilder.from_dataframe(df_display)
    
    # Selection Mode
    gb.configure_selection('single', use_checkbox=False)
    
    # Styling und Formatierung
    gb.configure_column("Smiles", hide=True) # Smiles verstecken wir, brauchen sie aber im Hintergrund
    gb.configure_column("Affinity_Target", header_name="Target (kcal)", type=["numericColumn"], precision=2)
    
    if has_reference:
        gb.configure_column("Affinity_Ref", header_name="Ref (kcal)", type=["numericColumn"], precision=2)
        gb.configure_column("Delta_Affinity", header_name="Specificity", type=["numericColumn"], precision=2, sort="asc")
        
        # Färbe die Delta Spalte (optionales JS, hier einfach gehalten)
        gb.configure_column("Delta_Rank", header_name="Rank Gain")

    gb.configure_grid_options(domLayout='normal')
    gridOptions = gb.build()

    grid_response = AgGrid(
        df_display,
        gridOptions=gridOptions,
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        fit_columns_on_grid_load=True,
        height=500,
        theme='streamlit'
    )

with col_viz:
    st.write("### Structure Preview")
    
    selected = grid_response['selected_rows']
    
    # Logik um die Auswahl sicher zu extrahieren
    if selected is not None and len(selected) > 0:
        # AgGrid gibt manchmal DataFrame, manchmal Liste von Dicts zurück
        if isinstance(selected, pd.DataFrame):
            selected_row = selected.iloc[0]
        else:
            selected_row = pd.DataFrame(selected).iloc[0]
            
        smiles = selected_row.get("Smiles", "")
        ligand_name = selected_row.get("Ligand", "Unknown")
        
        if smiles and smiles != "unknown":
            # RDKit Zeichnung
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(350, 350))
                    st.image(img, caption=f"Ligand: {ligand_name}")
                    
                    # Details anzeigen
                    st.info(f"**SMILES:** `{smiles}`")
                    if has_reference:
                        st.metric("Specificity Score", f"{selected_row['Delta_Affinity']:.2f} kcal/mol")
                else:
                    st.error("Invalid SMILES structure.")
            except Exception as e:
                st.error(f"Error rendering structure: {e}")
        else:
            st.warning("No SMILES string available for this ligand.")
    else:
        st.info("Select a ligand from the table to visualize it here.")
        # Platzhalter-Bild (optional)
        # st.image("https://placehold.co/300x300?text=Select+Ligand", use_column_width=True)