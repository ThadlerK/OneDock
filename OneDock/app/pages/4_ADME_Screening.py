import streamlit as st
import pandas as pd
import os
import time
import requests
from io import StringIO

st.title("ADME Screening")
st.markdown("")

# Ensure output directory exists
os.makedirs("data/results/output", exist_ok=True)

# Check if docking results exist
result_file = "data/results/docking_report_target.csv"

if not os.path.exists(result_file):
    st.warning("No docking results found. Please run the docking pipeline first.")
    if st.button("Go to Docking"):
        st.switch_page("pages/2_Docking.py")
    st.stop()

# Load docking results
df = pd.read_csv(result_file)

# Rename columns for clarity
if 'Ligand' in df.columns:
    df['Ligand_ID'] = df['Ligand']
if 'Smiles' in df.columns:
    df['SMILES'] = df['Smiles']

# Load reference results if available for specificity calculation
ref_file = "data/results/docking_report_reference.csv"
if os.path.exists(ref_file):
    df_ref = pd.read_csv(ref_file)
    if 'Ligand' in df_ref.columns:
        df_ref = df_ref.rename(columns={'Ligand': 'Ligand_ID'})
    
    df = df.merge(
        df_ref[['Ligand_ID', 'Affinity_kcal_mol']],
        on='Ligand_ID',
        how='left',
        suffixes=('', '_ref')
    )
    df['Specificity'] = df['Affinity_kcal_mol'] - df['Affinity_kcal_mol_ref']
    has_specificity = True
else:
    has_specificity = False

df_sorted = df.sort_values(by="Affinity_kcal_mol")

st.subheader("Select Ligands for ADME Screening")

# --- INITIALIZE SESSION STATE ---
if 'affinity_cutoff' not in st.session_state:
    st.session_state.affinity_cutoff = -3.0
if 'specificity_cutoff' not in st.session_state:
    st.session_state.specificity_cutoff = 0.0

# --- CUSTOM CUTOFF FILTERS ---
col1, col2 = st.columns(2)

with col1:
    affinity_cutoff = st.number_input(
        "Target Affinity ≤ (kcal/mol):",
        value=st.session_state.affinity_cutoff,
        step=0.5,
        key='adme_affinity'
    )
    st.session_state.affinity_cutoff = affinity_cutoff

with col2:
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Specificity ≤:",
            value=st.session_state.specificity_cutoff,
            step=0.5,
            key='adme_specificity'
        )
        st.session_state.specificity_cutoff = specificity_cutoff
    else:
        st.info("No reference data - specificity filter not available")
        specificity_cutoff = None

# Apply filters
if has_specificity and specificity_cutoff is not None:
    selected_df = df_sorted[
        (df_sorted['Affinity_kcal_mol'] <= affinity_cutoff) &
        (df_sorted['Specificity'] <= specificity_cutoff)
    ][['Ligand_ID', 'Affinity_kcal_mol', 'Specificity', 'SMILES']].copy()
else:
    selected_df = df_sorted[
        df_sorted['Affinity_kcal_mol'] <= affinity_cutoff
    ][['Ligand_ID', 'Affinity_kcal_mol', 'SMILES']].copy()

selected_ligands = selected_df['Ligand_ID'].tolist()

st.dataframe(selected_df, height=250)
st.info(f"Found {len(selected_ligands)} ligands matching criteria")

# --- SWISSADME SUBMISSION ---
st.subheader("Submit to SwissADME")

if selected_ligands:
    # Extract SMILES for submission (already in dataframe)
    smiles_data = []
    
    for ligand_id in selected_ligands:
        # Get SMILES from dataframe
        ligand_row = df_sorted[df_sorted['Ligand_ID'] == ligand_id]
        if not ligand_row.empty:
            smiles = ligand_row['SMILES'].values[0]
            smiles_data.append({
                'Ligand_ID': ligand_id,
                'SMILES': smiles
            })
    
    if smiles_data:
        smiles_df = pd.DataFrame(smiles_data)
        
        # Prepare SMILES list for SwissADME
        smiles_list = "\n".join([f"{row['SMILES']} {row['Ligand_ID']}" for _, row in smiles_df.iterrows()])
        
        st.text_area("SMILES for SwissADME (Copy this):", value=smiles_list, height=200)
        
        # Manual submission instructions
        st.markdown("""
        **How to submit to SwissADME:**
        1. Copy the SMILES list above
        2. Open [SwissADME](https://www.swissadme.ch/index.php) in a new tab
        3. Paste the SMILES into the input field
        4. Click "Run" and wait for results
        5. Download the results CSV
        6. Upload it below for integration with docking results
        """)
        
        if st.button("Open SwissADME in Browser", key="open_swissadme"):
            st.markdown("[Click here to open SwissADME](https://www.swissadme.ch/index.php)", unsafe_allow_html=True)
        
        # --- UPLOAD SWISSADME RESULTS ---
        st.subheader("Upload SwissADME Results")
        
        uploaded_adme = st.file_uploader(
            "Upload SwissADME results (.csv)",
            type=['csv']
        )
        
        if uploaded_adme:
            try:
                adme_df = pd.read_csv(uploaded_adme)

                # Save ADME results
                adme_output_path = "data/results/output/swissadme_results.csv"
                adme_df.to_csv(adme_output_path, index=False)
                
                # --- MERGE WITH DOCKING RESULTS ---

                # Merge dataframes (assuming 'Name' or 'Molecule Name' column in SwissADME matches Ligand_ID)
                # You may need to adjust the merge key depending on SwissADME output format
                if 'Name' in adme_df.columns:
                    merge_key = 'Name'
                elif 'Molecule Name' in adme_df.columns:
                    merge_key = 'Molecule Name'
                else:
                    merge_key = adme_df.columns[0]
                
                # Merge with docking results
                merged_df = df_sorted.merge(
                    adme_df,
                    left_on='Ligand_ID',
                    right_on=merge_key,
                    how='inner'
                )
                
                if not merged_df.empty:
                    
                    
                    # Mapping: Anzeige-Name -> tatsächlicher Spaltenname in SwissADME / Docking
                    column_mapping = {
                        'Ligand_ID': 'Ligand_ID',             # aus Docking
                        'SMILES': 'SMILES',                   # SMILES der Liganden
                        'Affinity_kcal_mol': 'Affinity_kcal_mol',
                        'Formula': 'Formula',
                        'Molecular Weight': 'MW',              # SwissADME: MW
                        'Num. H-bond acceptors': '#H-bond acceptors',
                        'Num. H-bond donors': '#H-bond donors',
                        'Consensus Log Po/w': 'Consensus Log P',
                        'GI absorption': 'GI absorption',
                        'BBB permeant': 'BBB permeant',
                        'Lipinski #violations': 'Lipinski #violations',
                        'Bioavailability Score': 'Bioavailability Score',
                        'Synthetic accessibility': 'Synthetic Accessibility'
                    }

                    display_df = {}
                    missing_columns = []

                    for display_name, source_name in column_mapping.items():
                        if source_name in merged_df.columns:
                            display_df[display_name] = merged_df[source_name]
                        else:
                            missing_columns.append(display_name)

                    if display_df:
                        display_df = pd.DataFrame(display_df)

                        # Nach Affinität sortieren (aufsteigend = beste zuerst)
                        if 'Affinity_kcal_mol' in display_df.columns:
                            display_df = display_df.sort_values('Affinity_kcal_mol')

                        # Tabelle mit begrenzter Höhe anzeigen (Scroll-Bar für weitere Liganden)
                        st.dataframe(display_df, use_container_width=True, height=250)

                        # Warnung für wirklich fehlende Felder (falls SwissADME-Export reduziert ist)
                        if missing_columns:
                            st.warning(
                                "Some requested columns are not present in the SwissADME file: "
                                + ", ".join(missing_columns)
                            )
                    else:
                        st.warning("Could not find expected columns in SwissADME results")
                        st.dataframe(merged_df)
                    
                    # Save combined results
                    combined_output = "data/results/output/combined_docking_adme.csv"
                    merged_df.to_csv(combined_output, index=False)

                    # Download button
                    csv = merged_df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="Download Combined Results",
                        data=csv,
                        file_name="combined_docking_adme.csv",
                        mime="text/csv"
                    )
                    
                    # --- FILTERING OPTIONS ---
                    st.subheader("Druglikeness")

                    # Kurzinfo zur Lipinski-Rule-of-Five
                    st.info(
                        "**Lipinski's Rule of Five** defines criteria for oral bioavailability:\n\n"
                        "- MW ≤ 500 Da\n"
                        "- H-bond donors ≤ 5\n"
                        "- H-bond acceptors ≤ 10\n"
                        "- cLogP ≤ 5\n\n"
                        "Compounds violating more than one criterion typically exhibit poor oral absorption, though exceptions exist."
                    )

                    # Slider fr die einzelnen Lipinski-Parameter
                    if all(col in merged_df.columns for col in ['MW', '#H-bond donors', '#H-bond acceptors', 'Consensus Log P']):
                        col1, col2 = st.columns(2)

                        with col1:
                            max_hbd = st.slider(
                                "Max. H-bond donors",
                                min_value=0,
                                max_value=15,
                                value=5,
                                step=1,
                                help="Lipinski threshold: ≤ 5"
                            )
                            max_hba = st.slider(
                                "Max. H-bond acceptors",
                                min_value=0,
                                max_value=20,
                                value=10,
                                step=1,
                                help="Lipinski threshold: ≤ 10"
                            )

                        with col2:
                            max_mw = st.slider(
                                "Max. molecular weight (Da)",
                                min_value=0,
                                max_value=1000,
                                value=500,
                                step=10,
                                help="Lipinski threshold: ≤ 500 Da"
                            )
                            max_logp = st.slider(
                                "Max. LogP",
                                min_value=-10.0,
                                max_value=10.0,
                                value=5.0,
                                step=0.1,
                                help="Lipinski threshold: ≤ 5"
                            )

                        mask = (
                            (merged_df['#H-bond donors'] <= max_hbd) &
                            (merged_df['#H-bond acceptors'] <= max_hba) &
                            (merged_df['MW'] <= max_mw) &
                            (merged_df['Consensus Log P'] <= max_logp)
                        )

                        filtered_df = merged_df[mask]

                        st.write(f"Filtered: {len(filtered_df)} compounds")

                        # sinnvolle Standardspalten fr die Ansicht
                        cols_to_show = [
                            'Ligand_ID',
                            'SMILES',
                            'Affinity_kcal_mol',
                            'MW',
                            '#H-bond donors',
                            '#H-bond acceptors',
                            'Consensus Log P',
                            'Lipinski #violations'
                        ]
                        existing_cols = [c for c in cols_to_show if c in filtered_df.columns]
                        st.dataframe(filtered_df[existing_cols], use_container_width=True, height=250)
                    else:
                        st.warning(
                            "Lipinski-Filter kann nicht berechnet werden, da erforderliche Spalten "
                            "(MW, #H-bond donors, #H-bond acceptors, Consensus Log P) im SwissADME-Export fehlen."
                        )
                    
                    # --- VISUALIZATION ---
                    st.subheader("Molecular Structures of Lead Compounds")
                    
                    # Molecular structure viewer
                    st.markdown("")
                    compound_to_view = st.selectbox(
                        "Select compound:",
                        options=filtered_df['Ligand_ID'].tolist() if 'Ligand_ID' in filtered_df.columns else [],
                        key="structure_viewer"
                    )
                    
                    if compound_to_view:
                        # Get SMILES for selected compound
                        compound_row = filtered_df[filtered_df['Ligand_ID'] == compound_to_view]
                        
                        if not compound_row.empty and 'SMILES' in compound_row.columns:
                            try:
                                from rdkit import Chem
                                from rdkit.Chem import Draw
                                from io import BytesIO
                                
                                smiles_str = compound_row['SMILES'].values[0]
                                mol = Chem.MolFromSmiles(smiles_str)
                                
                                if mol:
                                    # Generate 2D structure image
                                    img = Draw.MolToImage(mol, size=(400, 400))
                                    
                                    # Display structure
                                    col1, col2 = st.columns([1, 2])
                                    with col1:
                                        st.image(img, caption=f"{compound_to_view}", use_container_width=True)
                                    
                                    with col2:
                                        # Show compound properties
                                        compound_data = merged_df[merged_df['Ligand_ID'] == compound_to_view].iloc[0]
                                        st.markdown(f"**Compound:** {compound_to_view}")
                                        st.markdown(f"**SMILES:** `{smiles_str}`")
                                        st.markdown(f"**Binding Affinity:** {compound_data['Affinity_kcal_mol']:.2f} kcal/mol")
                                        
                                        if 'Molecular Weight' in compound_data:
                                            st.markdown(f"**Molecular Weight:** {compound_data['Molecular Weight']:.1f} Da")
                                        if 'LogP' in compound_data:
                                            st.markdown(f"**LogP:** {compound_data['LogP']:.2f}")
                                        if 'Lipinski #violations' in compound_data:
                                            st.markdown(f"**Lipinski Violations:** {int(compound_data['Lipinski #violations'])}")
                                else:
                                    st.error("Could not parse SMILES structure")
                            except ImportError:
                                st.warning("RDKit not installed. Install with: `pip install rdkit`")
                            except Exception as e:
                                st.error(f"Error displaying structure: {str(e)}")
                    
                    # Create scatter plot if relevant columns exist
                    if 'Molecular Weight' in merged_df.columns and 'TPSA' in merged_df.columns:
                        import plotly.express as px
                        
                        fig = px.scatter(
                            merged_df,
                            x='Molecular Weight',
                            y='TPSA',
                            color='Affinity_kcal_mol',
                            hover_data=['Ligand_ID'],
                            title='Molecular Weight vs TPSA (colored by binding affinity)',
                            labels={'TPSA': 'Topological Polar Surface Area (Ų)'}
                        )
                        st.plotly_chart(fig, use_container_width=True)
                    
                else:
                    st.warning("Could not merge results. Check that ligand names match between files.")
                    st.write("Docking ligand names:", df_sorted['Ligand_ID'].head().tolist())
                    st.write("ADME column names:", adme_df.columns.tolist())
                
            except Exception as e:
                st.error(f"Error processing SwissADME results: {str(e)}")
                st.info("Make sure you uploaded the correct CSV file from SwissADME.")

else:
    st.warning("Please select ligands to proceed with ADME screening.")

# --- NAVIGATION ---
col1, col2 = st.columns(2)
with col1:
    if st.button("← Back to Results"):
        st.switch_page("pages/3_Results.py")
