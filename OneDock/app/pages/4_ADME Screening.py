import streamlit as st
import pandas as pd
import os
import time
import requests
from io import StringIO
from utils import load_config

st.title("ADME Screening")
st.markdown("### SwissADME-based Property Prediction")

# Display information about SwissADME
st.info("""
**SwissADME** is a web-based tool that uses a variety of predictive models to compute physicochemical descriptors of small molecules from their structure.

Reference: Daina, A., Michielin, O., & Zoete, V. (2017). SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. *Scientific reports, 7*(1), 42717.
""")

st.set_page_config(page_title = "ADME Screening")
# Ensure output directory exists
os.makedirs("data/results/output", exist_ok=True)

# Check if docking results exist
config = load_config()

lib_name = config.get("library_name", '')
result_file = f"data/results/docking_report_target_{lib_name}.csv"

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
ref_file = f"data/results/docking_report_reference_{lib_name}.csv"
if os.path.exists(ref_file):
    df_ref = pd.read_csv(ref_file)
    if 'Ligand' in df_ref.columns:
        df_ref = df_ref.rename(columns={'Ligand': 'Ligand_ID'})
    
    # Sort by affinity to assign ranks
    df_target_sorted = df.sort_values('Affinity_kcal_mol').reset_index(drop=True)
    df_target_sorted['Rank_Target'] = df_target_sorted.index + 1
    
    df_ref_sorted = df_ref.sort_values('Affinity_kcal_mol').reset_index(drop=True)
    df_ref_sorted['Rank_Ref'] = df_ref_sorted.index + 1
    
    # Merge ranks and reference affinity
    df = df.merge(df_target_sorted[['Ligand_ID', 'Rank_Target']], on='Ligand_ID', how='left')
    df = df.merge(df_ref_sorted[['Ligand_ID', 'Rank_Ref']], on='Ligand_ID', how='left')
    df = df.merge(
        df_ref[['Ligand_ID', 'Affinity_kcal_mol']],
        on='Ligand_ID',
        how='left',
        suffixes=('', '_ref')
    )
    df['Specificity'] = df['Affinity_kcal_mol'] - df['Affinity_kcal_mol_ref']
    df['Rank_Gain'] = df['Rank_Ref'] - df['Rank_Target']
    has_specificity = True
    has_rank_gain = True
else:
    has_specificity = False
    has_rank_gain = False

df_sorted = df.sort_values(by="Affinity_kcal_mol")

st.subheader("Select Ligands for ADME Screening")

# --- INITIALIZE SESSION STATE ---
if 'affinity_cutoff' not in st.session_state:
    st.session_state.affinity_cutoff = -3.0
if 'specificity_cutoff' not in st.session_state:
    st.session_state.specificity_cutoff = 0.0
if 'adme_merged_df' not in st.session_state:
    st.session_state.adme_merged_df = None
if 'adme_filter_values' not in st.session_state:
    st.session_state.adme_filter_values = {
        'max_hbd': 5,
        'max_hba': 10,
        'max_mw': 500,
        'max_logp': 5.0
    }

# --- CUSTOM CUTOFF FILTERS ---
col1, col2, col3 = st.columns(3)

with col1:
    affinity_cutoff = st.number_input(
        "Max Target Affinity (kcal/mol):",
        value=st.session_state.affinity_cutoff,
        step=0.5,
        key='adme_affinity',
        help="Show only ligands with affinity ≤ this value."
    )
    st.session_state.affinity_cutoff = affinity_cutoff

with col2:
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Max Specificity:",
            value=st.session_state.specificity_cutoff,
            step=0.5,
            key='adme_specificity',
            help="Show only ligands with specificity ≤ this value."
        )
        st.session_state.specificity_cutoff = specificity_cutoff
    else:
        st.warning("⚠️ No reference data - specificity filter not available")
        specificity_cutoff = None

with col3:
    if has_rank_gain:
        min_rank_gain = int(df['Rank_Gain'].min())
        max_rank_gain = int(df['Rank_Gain'].max())
        rank_gain_cutoff = st.number_input(
            "Min Rank Gain:",
            min_value=min_rank_gain,
            max_value=max_rank_gain,
            value=st.session_state.get('rank_gain_cutoff', min_rank_gain),
            step=1,
            key='adme_rank_gain',
            help="Show only ligands with rank gain ≥ this value."
        )
        st.session_state.rank_gain_cutoff = rank_gain_cutoff
    else:
        st.warning("⚠️ No reference data - rank gain filter not available")
        rank_gain_cutoff = None

# Apply filters
selected_df = df_sorted[df_sorted['Affinity_kcal_mol'] <= affinity_cutoff].copy()

if has_specificity and specificity_cutoff is not None:
    selected_df = selected_df[selected_df['Specificity'] <= specificity_cutoff]

if has_rank_gain and rank_gain_cutoff is not None:
    selected_df = selected_df[selected_df['Rank_Gain'] >= rank_gain_cutoff]

# Select display columns
display_cols = ['Ligand_ID', 'Affinity_kcal_mol']
if has_specificity:
    display_cols.append('Specificity')
if has_rank_gain:
    display_cols.append('Rank_Gain')
display_cols.append('SMILES')

selected_df = selected_df[display_cols].copy()

selected_ligands = selected_df['Ligand_ID'].tolist()

st.dataframe(selected_df, height=250)

# --- CHECK IF ADME DATA IS ALREADY LOADED ---
if st.session_state.adme_merged_df is not None:
    
    col_clear1, col_clear2 = st.columns([3, 1])
    with col_clear2:
        if st.button("Remove ADME Data & Upload New", type="secondary"):
            st.session_state.adme_merged_df = None
            st.session_state.adme_filter_values = {
                'max_hbd': 5,
                'max_hba': 10,
                'max_mw': 500,
                'max_logp': 5.0
            }
            st.rerun()
    
    merged_df = st.session_state.adme_merged_df
    
    # Display merged data
    st.subheader("ADME Results")
    
    # Mapping for display
    column_mapping = {
        'Ligand_ID': 'Ligand_ID',
        'SMILES': 'SMILES',
        'Affinity_kcal_mol': 'Affinity_kcal_mol',
        'Formula': 'Formula',
        'Molecular Weight': 'MW',
        'Num. H-bond acceptors': '#H-bond acceptors',
        'Num. H-bond donors': '#H-bond donors',
        'Consensus Log Po/w': 'Consensus Log P',
        'TPSA': 'TPSA',
        'Rotatable Bonds': '#Rotatable bonds',
        'Aromatic Rings': '#Aromatic heavy atoms',
        'GI absorption': 'GI absorption',
        'BBB permeant': 'BBB permeant',
        'Lipinski #violations': 'Lipinski #violations',
        'PAINS alerts': 'PAINS #alerts',
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
        if 'Affinity_kcal_mol' in display_df.columns:
            display_df = display_df.sort_values('Affinity_kcal_mol')
        st.dataframe(display_df, use_container_width=True, height=250)
        
        if missing_columns:
            st.warning(f"Some columns are not present: {', '.join(missing_columns)}")
    
    # Parameter descriptions
    with st.expander("Parameter Information"):
        st.markdown("""
        **Molecular Weight (MW):** Mass of the molecule in Daltons. Drug-like compounds typically have MW&nbsp;≤&nbsp;500&nbsp;Da.
        
        **H-bond Donors:** Number of hydrogen atoms attached to nitrogen or oxygen atoms that can donate hydrogen bonds. Lipinski:&nbsp;≤&nbsp;5.
        
        **H-bond Acceptors:** Number of nitrogen or oxygen atoms that can accept hydrogen bonds. Lipinski:&nbsp;≤&nbsp;10.
        
        **Consensus Log P:** Average of five different calculation methods for the octanol-water partition coefficient. Measures lipophilicity. Lipinski:&nbsp;≤&nbsp;5.
        
        **TPSA (Topological Polar Surface Area):** Sum of surfaces of polar atoms (usually oxygens and nitrogens). Influences membrane permeability. Good oral bioavailability:&nbsp;<&nbsp;140&nbsp;Ų.
        
        **Rotatable Bonds:** Number of bonds that allow free rotation. Too many reduce oral bioavailability. Recommended:&nbsp;<&nbsp;10.
        
        **Aromatic Rings:** Number of aromatic heavy atoms. Important for protein binding but too many can reduce bioavailability. Optimal:&nbsp;1-5&nbsp;rings.
        
        **GI Absorption:** Predicted gastrointestinal absorption (High/Low).
        
        **BBB Permeant:** Blood-Brain Barrier permeability (Yes/No).
        
        **Lipinski #violations:** Number of Lipinski Rule of Five violations. 0-1&nbsp;violations indicate drug-like properties.
        
        **PAINS alerts:** Pan-Assay Interference Compounds - structural patterns that frequently cause false-positive results in assays. Should be&nbsp;0 for reliable compounds.
        
        **Bioavailability Score:** Probability of good oral bioavailability. SwissADME assigns discrete values: 0.55&nbsp;or&nbsp;0.56&nbsp;=&nbsp;Excellent, 0.17&nbsp;=&nbsp;Medium, 0.11&nbsp;=&nbsp;Low, <&nbsp;0.11&nbsp;=&nbsp;Poor. Higher is better.
        
        **Synthetic Accessibility:** Estimated ease of synthesis (1&nbsp;=&nbsp;very easy, 10&nbsp;=&nbsp;very difficult).
        """, unsafe_allow_html=True)
    
    # Download button
    csv = merged_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Combined Results",
        data=csv,
        file_name="combined_docking_adme.csv",
        mime="text/csv"
    )
    
    # --- DRUGLIKENESS FILTERING ---
    st.markdown("---")
    st.subheader("Druglikeness Filtering")
    
    with st.container(border=True):
        st.markdown(
            "**Lipinski's Rule of Five** defines criteria for oral bioavailability:\n\n"
            "- MW ≤ 500 Da\n"
            "- H-bond donors ≤ 5\n"
            "- H-bond acceptors ≤ 10\n"
            "- cLogP ≤ 5\n\n"
            "Compounds violating more than one criterion typically exhibit poor oral absorption."
        )
    
    if all(col in merged_df.columns for col in ['MW', '#H-bond donors', '#H-bond acceptors', 'Consensus Log P']):
        col1, col2 = st.columns(2)
        
        with col1:
            max_hbd = st.slider(
                "Max. H-bond donors",
                min_value=0,
                max_value=15,
                value=st.session_state.adme_filter_values['max_hbd'],
                step=1,
                help="Lipinski threshold: ≤ 5",
                key="slider_hbd"
            )
            st.session_state.adme_filter_values['max_hbd'] = max_hbd
            
            max_hba = st.slider(
                "Max. H-bond acceptors",
                min_value=0,
                max_value=20,
                value=st.session_state.adme_filter_values['max_hba'],
                step=1,
                help="Lipinski threshold: ≤ 10",
                key="slider_hba"
            )
            st.session_state.adme_filter_values['max_hba'] = max_hba
        
        with col2:
            max_mw = st.slider(
                "Max. molecular weight (Da)",
                min_value=0,
                max_value=1000,
                value=st.session_state.adme_filter_values['max_mw'],
                step=10,
                help="Lipinski threshold: ≤ 500 Da",
                key="slider_mw"
            )
            st.session_state.adme_filter_values['max_mw'] = max_mw
            
            max_logp = st.slider(
                "Max. LogP",
                min_value=-10.0,
                max_value=10.0,
                value=st.session_state.adme_filter_values['max_logp'],
                step=0.1,
                help="Lipinski threshold: ≤ 5",
                key="slider_logp"
            )
            st.session_state.adme_filter_values['max_logp'] = max_logp
        
        mask = (
            (merged_df['#H-bond donors'] <= max_hbd) &
            (merged_df['#H-bond acceptors'] <= max_hba) &
            (merged_df['MW'] <= max_mw) &
            (merged_df['Consensus Log P'] <= max_logp)
        )
        
        filtered_df = merged_df[mask]
        
        st.write(f"**Filtered: {len(filtered_df)} compounds match the criteria**")
        
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
        
        # --- VISUALIZATION ---
        st.markdown("---")
        st.subheader("Molecular Structures of Lead Compounds")
        
        compound_to_view = st.selectbox(
            "Select compound:",
            options=filtered_df['Ligand_ID'].tolist() if 'Ligand_ID' in filtered_df.columns else [],
            key="structure_viewer"
        )
        
        if compound_to_view:
            compound_row = filtered_df[filtered_df['Ligand_ID'] == compound_to_view]
            if not compound_row.empty and 'SMILES' in compound_row.columns:
                smiles_str = compound_row['SMILES'].values[0]
                
                try:
                    from rdkit import Chem
                    from rdkit.Chem import Draw
                    
                    mol = Chem.MolFromSmiles(smiles_str)
                    if mol:
                        img = Draw.MolToImage(mol, size=(400, 300))
                        st.image(img, caption=f"Structure of {compound_to_view}")
                        
                        # Display properties
                        st.markdown("**Properties:**")
                        prop_cols = st.columns(3)
                        with prop_cols[0]:
                            if 'MW' in compound_row.columns:
                                st.metric("MW", f"{compound_row['MW'].values[0]:.2f} Da")
                        with prop_cols[1]:
                            if 'Consensus Log P' in compound_row.columns:
                                st.metric("LogP", f"{compound_row['Consensus Log P'].values[0]:.2f}")
                        with prop_cols[2]:
                            if 'Affinity_kcal_mol' in compound_row.columns:
                                st.metric("Affinity", f"{compound_row['Affinity_kcal_mol'].values[0]:.2f} kcal/mol")
                    else:
                        st.error("Could not parse SMILES")
                except Exception as e:
                    st.error(f"Error generating structure: {e}")
    else:
        st.warning("Required columns for Lipinski filtering are missing in the SwissADME data.")
    
    st.stop()

# --- SWISSADME SUBMISSION (Only shown if no data loaded) ---
st.markdown("")
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
                    
                    # Store in session state
                    st.session_state.adme_merged_df = merged_df
                    
                    st.success("ADME data successfully loaded and merged")
                    st.rerun()
                
                else:
                    st.error("No matches found between docking results and SwissADME data.")
            
            except Exception as e:
                st.error(f"Error reading SwissADME file: {e}")
else:
    st.write("No ligands selected. Adjust filters to select ligands for ADME screening.")
