import streamlit as st
import pandas as pd
import numpy as np
import glob
import os
import requests
import urllib.parse
import time
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
from rdkit import Chem
from rdkit.Chem import Draw
from utils import load_config

st.set_page_config(layout="wide", page_title="Summary Report")

st.title("Summary Report")

# Add reload button
if st.button("Reload Data"):
    st.cache_data.clear()
    st.rerun()

# --- HELPER FUNCTIONS ---
@st.cache_data
def fetch_pubchem_data(smiles):
    """
    Queries PubChem ONCE to get Common Name, IUPAC, Description, and Link.
    Returns a dictionary or None.
    """
    if not smiles or len(str(smiles)) < 2: 
        return None

    data = {
        "CommonName": None,
        "IUPAC": None,
        "Description": None,
        "Link": None
    }

    try:
        encoded_smiles = urllib.parse.quote(smiles)
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        
        # 1. Get CID (Compound ID)
        url_cid = f"{base_url}/compound/smiles/{encoded_smiles}/cids/JSON"
        resp = requests.get(url_cid, timeout=2)
        if resp.status_code != 200:
            return None
        
        cid = resp.json()['IdentifierList']['CID'][0]
        data["Link"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

        # 2. Get Properties (Title, IUPAC Name) in ONE call
        url_props = f"{base_url}/compound/cid/{cid}/property/Title,IUPACName/JSON"
        resp_props = requests.get(url_props, timeout=2)
        if resp_props.status_code == 200:
            props = resp_props.json()['PropertyTable']['Properties'][0]
            data["CommonName"] = props.get('Title', f"CID: {cid}")
            data["IUPAC"] = props.get('IUPACName')

        # 3. Get Description (Separate call usually required)
        url_desc = f"{base_url}/compound/cid/{cid}/description/JSON"
        resp_desc = requests.get(url_desc, timeout=2)
        if resp_desc.status_code == 200:
            desc_data = resp_desc.json()
            if 'InformationList' in desc_data and 'Information' in desc_data['InformationList']:
                # Get the first available description
                data["Description"] = desc_data['InformationList']['Information'][0].get('Description')

        return data
            
    except Exception:
        return None

def get_pocket_residues_from_pose(receptor_pdbqt, ligand_pdbqt, cutoff=3.5):
    """Extract residues within cutoff distance from docked ligand"""
    try:
        # Read receptor coordinates
        receptor_atoms = []
        with open(receptor_pdbqt, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        chain = line[21].strip()
                        resnum = int(line[22:26].strip())
                        resname = line[17:20].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        receptor_atoms.append({
                            'chain': chain if chain else 'A',
                            'resnum': resnum,
                            'resname': resname,
                            'coords': np.array([x, y, z])
                        })
                    except (ValueError, IndexError):
                        continue
        
        # Read ligand coordinates
        ligand_coords = []
        with open(ligand_pdbqt, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append(np.array([x, y, z]))
                    except (ValueError, IndexError):
                        continue
        
        if not ligand_coords or not receptor_atoms:
            return "N/A"
        
        # Find residues within cutoff
        nearby_residues = set()
        for lig_coord in ligand_coords:
            for rec_atom in receptor_atoms:
                distance = np.linalg.norm(lig_coord - rec_atom['coords'])
                if distance <= cutoff:
                    res_id = f"{rec_atom['chain']}:{rec_atom['resname']}{rec_atom['resnum']}"
                    nearby_residues.add(res_id)
        
        if nearby_residues:
            return ", ".join(sorted(nearby_residues))
        else:
            return "N/A"
            
    except Exception as e:
        return f"Error: {str(e)}"

# --- LOAD DATA ---
config = load_config()

lib_name = config.get("library_name", '')
RUN_NAME = config.get("run_name", 'default_run')
TARGET_FILE = f"data/results/docking_report_target_{lib_name}_{RUN_NAME}.csv"
REF_FILE = f"data/results/docking_report_reference_{lib_name}_{RUN_NAME}.csv"


@st.cache_data
def load_summary_data(target_path, ref_path):
    if not os.path.exists(target_path):
        return None
    
    # Load target docking results
    df = pd.read_csv(target_path)
    df['Affinity_kcal_mol'] = pd.to_numeric(df['Affinity_kcal_mol'], errors='coerce')

    df['Rank_Target'] = df['Affinity_kcal_mol'].rank(method='min', ascending=True)
    
    # Rename columns
    df = df.rename(columns={
        'Affinity_kcal_mol': 'Affinity_Target',
        'Ligand': 'Ligand'
    })
    
    # Round affinity to 2 decimal places
    df['Affinity_Target'] = df['Affinity_Target'].round(2)
    
    # Extract rank from ligand name if present (e.g., lig_00000_rank1)
    if 'Ligand' in df.columns:
        df['Rank_String'] = df['Ligand'].str.extract(r'_rank(\d+)$')[0].fillna('1').astype(int)
    
    # Load reference if available and calculate specificity & rank gain
    if os.path.exists(ref_path):
        df_ref = pd.read_csv(ref_path)
    

        df_ref['Affinity_kcal_mol'] = pd.to_numeric(df_ref['Affinity_kcal_mol'], errors='coerce')
        df_ref['Rank_Ref'] = df_ref['Affinity_kcal_mol'].rank(method='min', ascending=True)
        

        df_ref = df_ref.rename(columns={
            'Affinity_kcal_mol': 'Affinity_Ref',
            'Ligand': 'Ligand'
        })
        
        ref_subset = df_ref[['Ligand', 'Affinity_Ref', 'Rank_Ref']]
        
        # merge
        df = df.merge(ref_subset, on='Ligand', how='left')
        

        df['Specificity'] = (df['Affinity_Target'] - df['Affinity_Ref']).round(2)
        df['Rank_Gain'] = df['Rank_Ref'] - df['Rank_Target']
    else:
        df['Specificity'] = None
        df['Rank_Gain'] = None

    cols_to_fix = ['Rank_Target', 'Rank_Ref', 'Rank_Gain']

    for col in cols_to_fix:
        if col in df.columns:
            df[col] = df[col].astype('Int64')
    
    # Load ADME data if available
    adme_file = "data/results/output/swissadme_results.csv"
    if os.path.exists(adme_file):
        adme_df = pd.read_csv(adme_file)
        # Determine merge key
        if 'Molecule' in adme_df.columns:
            merge_key = 'Molecule'
        elif 'Name' in adme_df.columns:
            merge_key = 'Name'
        else:
            merge_key = adme_df.columns[0]
        
        df = df.merge(adme_df, left_on='Ligand', right_on=merge_key, how='left')
    
    # # Calculate pocket residues for each ligand
    # receptor_pdbqt = "data/interim/target_prep.pdbqt"
    # pocket_residues_list = []
    
    # for idx, row in df.iterrows():
    #     ligand_id = row['Ligand']
    #     ligand_pdbqt = f"data/results/target/{RUN_NAME}/poses/{ligand_id}_docked.pdbqt"
        
    #     if os.path.exists(receptor_pdbqt) and os.path.exists(ligand_pdbqt):
    #         residues = get_pocket_residues_from_pose(receptor_pdbqt, ligand_pdbqt, cutoff=3.5)
    #         pocket_residues_list.append(residues)
    #     else:
    #         pocket_residues_list.append("N/A")
    
    # df['Pocket_Residues'] = pocket_residues_list
    
    # Load PoseBusters results if available
    pb_results_file = "data/results/posebusters/posebusters_validation.csv"
    if os.path.exists(pb_results_file):
        pb_df = pd.read_csv(pb_results_file)
        pb_df['PoseBusters_Score'] = (pb_df['quality_score'] * 100).round(1)
        df = df.merge(pb_df[['ligand_id', 'PoseBusters_Score']], left_on='Ligand', right_on='ligand_id', how='left')
        df = df.drop(columns=['ligand_id'], errors='ignore')
    
    # Calculate Lipinski Rule of Five compliance
    from rdkit.Chem import Descriptors
    ro5_status_list = []
    
    for idx, row in df.iterrows():
        if 'Smiles' in row and row['Smiles']:
            try:
                mol = Chem.MolFromSmiles(row['Smiles'])
                if mol:
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    hbd = Descriptors.NumHDonors(mol)
                    hba = Descriptors.NumHAcceptors(mol)
                    
                    violations = 0
                    if mw > 500: violations += 1
                    if logp > 5: violations += 1
                    if hbd > 5: violations += 1
                    if hba > 10: violations += 1
                    
                    if violations == 0:
                        ro5_status_list.append('Pass')
                    elif violations == 1:
                        ro5_status_list.append('Pass (1 violation)')
                    else:
                        ro5_status_list.append('Fail')
                else:
                    ro5_status_list.append('N/A')
            except:
                ro5_status_list.append('N/A')
        else:
            ro5_status_list.append('N/A')
    
    df['Lipinski_RO5'] = ro5_status_list

    # include mmpbsa results
    results_files = glob.glob("data/results/mmpbsa/*/FINAL_RESULTS_MMPBSA.dat")

    mmpbsa_data = []

    if results_files:
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
                            
                            mmpbsa_data.append({
                                "Ligand": lig_id, 
                                "Rank": rank_id, 
                                "Delta G (MMGBSA)": val
                            })
                            break
            except:
                pass
        
    if mmpbsa_data and df is not None:
        df_mmpbsa = pd.DataFrame(mmpbsa_data)
        # B. Handle multiple ranks (Optional but recommended)
        df_mmpbsa_best = df_mmpbsa.sort_values("Delta G (MMGBSA)").groupby("Ligand").first().reset_index()

        # C. Merge into main dataframe
        df = df.merge(df_mmpbsa_best[['Ligand', 'Delta G (MMGBSA)']], on='Ligand', how='left')
                    
    
    # Sort by affinity (best first)
    df = df.sort_values('Affinity_Target', ascending=True)
    
    return df

df = load_summary_data(TARGET_FILE, REF_FILE)

if df is None:
    st.warning("No docking results found. Please run the docking pipeline first.")
    if st.button("Go to Docking"):
        st.switch_page("pages/2_Docking.py")
    st.stop()

# --- INITIALIZE SESSION STATE ---
if 'affinity_cutoff' not in st.session_state:
    st.session_state.affinity_cutoff = 0.0
if 'specificity_cutoff' not in st.session_state:
    st.session_state.specificity_cutoff = 0.0
if 'rank_filter' not in st.session_state:
    st.session_state.rank_filter = 1

# Ensure lipinski_filter is always defined
if 'lipinski_filter' not in locals():
    lipinski_filter = "All"

# --- MAIN FILTER CONTROLS ---
col1, col2, col3 = st.columns(3)

with col1:
    affinity_cutoff = st.number_input(
        "Max Target Affinity (kcal/mol):",
        value=0.0,
        step=0.5,
        key='summary_affinity',
        help="Show only ligands with affinity ≤ this value."
    )

with col2:
    has_specificity = df['Specificity'].notna().any()
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Max Specificity:",
            value=0.0,
            step=0.5,
            key='summary_specificity',
            help="Show only ligands with specificity ≤ this value."
        )
    else:
        st.warning("⚠️ No reference data - specificity filter not available")
        specificity_cutoff = None

with col3:
    has_rank_gain = df['Rank_Gain'].notna().any()
    if has_rank_gain:
        min_rank_gain = int(df['Rank_Gain'].min())
        max_rank_gain = int(df['Rank_Gain'].max())
        rank_gain_cutoff = st.number_input(
            "Min Rank Gain:",
            min_value=min_rank_gain,
            max_value=max_rank_gain,
            value=min_rank_gain,
            step=1,
            key='summary_rank_gain',
            help="Show only ligands with rank gain ≥ this value."
        )
    else:
        st.warning("⚠️ No reference data - rank gain filter not available")
        rank_gain_cutoff = None

# --- ADVANCED ADME FILTERS ---

# --- ADVANCED ADME FILTERS (Dropdown) ---

# --- ADVANCED ADME FILTERS (Alle nebeneinander) ---
with st.expander("Advanced Filters"):
    adme_col1, adme_col2, adme_col3, adme_col4 = st.columns(4)
    with adme_col1:
        max_tpsa = st.number_input(
            "Max TPSA (Ų):",
            min_value=0.0,
            max_value=300.0,
            value=140.0,
            step=10.0,
            help="Topological Polar Surface Area. Recommended: < 140 Ų"
        )
        max_rotatable = st.number_input(
            "Max Rotatable Bonds:",
            min_value=0,
            max_value=20,
            value=10,
            step=1,
            help="Recommended: < 10"
        )
        max_aromatic = st.number_input(
            "Max Aromatic Atoms:",
            min_value=0,
            max_value=50,
            value=30,
            step=5,
            help="Optimal: 1-5 aromatic rings"
        )
    with adme_col2:
        max_pains = st.number_input(
            "Max PAINS Alerts:",
            min_value=0,
            max_value=10,
            value=0,
            step=1,
            help="Should be 0 for reliable compounds"
        )
        gi_filter = st.selectbox(
            "GI Absorption:",
            options=["All", "High", "Low"],
            help="Gastrointestinal absorption"
        )
        bbb_filter = st.selectbox(
            "BBB Permeant:",
            options=["All", "Yes", "No"],
            help="Blood-Brain Barrier permeability"
        )
    with adme_col3:
        min_bioavail = st.number_input(
            "Min Bioavailability Score:",
            min_value=0.0,
            max_value=1.0,
            value=0.11,
            step=0.01,
            help="0.55 = Good, 0.17 = Medium, 0.11 = Low"
        )
        max_synth_access = st.number_input(
            "Max Synthetic Accessibility:",
            min_value=1.0,
            max_value=10.0,
            value=10.0,
            step=0.1,
            help="Lower is easier to synthesize (1=easy, 10=difficult)"
        )
    with adme_col4:
        min_posebusters = st.number_input(
            "Min PoseBusters Quality Score (%):",
            min_value=0.0,
            max_value=100.0,
            value=0.0,
            step=1.0,
            help="80-100%: High quality pose"
        )
        lipinski_filter = st.selectbox(
            "Lipinski Rule of Five:",
            options=["All", "Pass", "Fail"],
            help="MW≤500, HBD≤5, HBA≤10, cLogP≤5. >1 violation: poor oral absorption."
        )

# --- APPLY FILTERS ---
# Start with base affinity filter
df = df[df['Affinity_Target'] <= affinity_cutoff]

# Apply specificity filter
if has_specificity and specificity_cutoff is not None:
    df = df[df['Specificity'] <= specificity_cutoff]

# Apply rank gain filter
if has_rank_gain and rank_gain_cutoff is not None:
    df = df[df['Rank_Gain'] >= rank_gain_cutoff]

# Apply ADME filters

def filter_with_na(df, col, op, value):
    if col not in df.columns:
        return df
    mask = df[col].isna() | op(df[col], value)
    return df[mask]

import operator
df = filter_with_na(df, 'TPSA', operator.le, max_tpsa)
df = filter_with_na(df, '#Rotatable bonds', operator.le, max_rotatable)
df = filter_with_na(df, '#Aromatic heavy atoms', operator.le, max_aromatic)
df = filter_with_na(df, 'PAINS #alerts', operator.le, max_pains)
if gi_filter != "All" and 'GI absorption' in df.columns:
    mask = df['GI absorption'].isna() | (df['GI absorption'] == gi_filter)
    df = df[mask]
if bbb_filter != "All" and 'BBB permeant' in df.columns:
    mask = df['BBB permeant'].isna() | (df['BBB permeant'] == bbb_filter)
    df = df[mask]
df = filter_with_na(df, 'Bioavailability Score', operator.ge, min_bioavail)
df = filter_with_na(df, 'Synthetic accessibility', operator.le, max_synth_access)

# Apply PoseBusters Score filter
if min_posebusters is not None:
    if 'PoseBusters (%)' in df.columns:
        mask = df['PoseBusters (%)'].isna() | (df['PoseBusters (%)'].astype(float) >= min_posebusters)
        df = df[mask]
    elif 'PoseBusters_Score' in df.columns:
        mask = df['PoseBusters_Score'].isna() | (df['PoseBusters_Score'].astype(float) >= min_posebusters)
        df = df[mask]

# Apply Lipinski filter
if (('Lipinski RO5' in df.columns) or ('Lipinski_RO5' in df.columns)) and lipinski_filter != "All":
    lip_col = 'Lipinski RO5' if 'Lipinski RO5' in df.columns else 'Lipinski_RO5'
    if lipinski_filter == "Pass":
        df = df[df[lip_col].str.startswith('Pass')]
    elif lipinski_filter == "Fail":
        df = df[df[lip_col] == 'Fail']

# --- DISPLAY SUMMARY TABLE ---
st.subheader("Summary Table")

# Prepare display dataframe - base columns
display_columns = ['Ligand', 'Affinity_Target', 'Specificity', 'Rank_Gain']

# Add optional columns
if 'PoseBusters_Score' in df.columns:
    display_columns.append('PoseBusters_Score')
if 'Lipinski_RO5' in df.columns:
    display_columns.append('Lipinski_RO5')

# Add ADME columns if available
adme_cols = ['TPSA', '#Rotatable bonds', '#Aromatic heavy atoms', 'PAINS #alerts', 
             'GI absorption', 'BBB permeant', 'Bioavailability Score', 'Synthetic accessibility']
for col in adme_cols:
    if col in df.columns:
        display_columns.append(col)

# Add SMILES at the end
if 'Smiles' in df.columns:
    display_columns.append('Smiles')
if 'Delta G (MMGBSA)' in df.columns:
    display_columns.append('Delta G (MMGBSA)')

df_display = df[display_columns].copy()

# Rename for better display
rename_dict = {
    'Affinity_Target': 'Target (kcal/mol)',
    'Pocket_Residues': 'Pocket Residues',
    'TPSA': 'TPSA (Ų)',
    '#Rotatable bonds': 'Rotatable Bonds',
    '#Aromatic heavy atoms': 'Aromatic Atoms',
    'PAINS #alerts': 'PAINS',
    'GI absorption': 'GI Absorption',
    'BBB permeant': 'BBB Permeant',
    'Bioavailability Score': 'Bioavailability',
    'Synthetic accessibility': 'Synth. Access.'
}
if 'PoseBusters_Score' in df_display.columns:
    rename_dict['PoseBusters_Score'] = 'PoseBusters (%)'
if 'Lipinski_RO5' in df_display.columns:
    rename_dict['Lipinski_RO5'] = 'Lipinski RO5'

df_display = df_display.rename(columns=rename_dict)

# Format all float columns with dot as decimal separator for display and CSV export
for col in df_display.select_dtypes(include=['float', 'float64']).columns:
    df_display[col] = df_display[col].map(lambda x: f"{x:.3f}" if pd.notnull(x) else "")

# Configure AgGrid
gb = GridOptionsBuilder.from_dataframe(df_display)
gb.configure_selection(selection_mode='single', use_checkbox=False)
gb.configure_column("Ligand", width=120, pinned='left')
gb.configure_column("Target (kcal/mol)", width=150, type=["numericColumn", "numberColumnFilter"])
gb.configure_column("Specificity", width=120, type=["numericColumn", "numberColumnFilter"])
if 'Rank Gain' in df_display.columns:
    gb.configure_column("Rank Gain", width=120, type=["numericColumn", "numberColumnFilter"])
# gb.configure_column("Pocket Residues", width=400, wrapText=True, autoHeight=True)

if 'PoseBusters (%)' in df_display.columns:
    gb.configure_column("PoseBusters (%)", width=130, type=["numericColumn", "numberColumnFilter"])
if 'Lipinski RO5' in df_display.columns:
    gb.configure_column("Lipinski RO5", width=120)

# Configure ADME columns
if 'TPSA (Ų)' in df_display.columns:
    gb.configure_column("TPSA (Ų)", width=110, type=["numericColumn", "numberColumnFilter"])
if 'Rotatable Bonds' in df_display.columns:
    gb.configure_column("Rotatable Bonds", width=130, type=["numericColumn", "numberColumnFilter"])
if 'Aromatic Atoms' in df_display.columns:
    gb.configure_column("Aromatic Atoms", width=130, type=["numericColumn", "numberColumnFilter"])
if 'PAINS' in df_display.columns:
    gb.configure_column("PAINS", width=100, type=["numericColumn", "numberColumnFilter"])
if 'GI Absorption' in df_display.columns:
    gb.configure_column("GI Absorption", width=120)
if 'BBB Permeant' in df_display.columns:
    gb.configure_column("BBB Permeant", width=120)
if 'Bioavailability' in df_display.columns:
    gb.configure_column("Bioavailability", width=130, type=["numericColumn", "numberColumnFilter"])
if 'Synth. Access.' in df_display.columns:
    gb.configure_column("Synth. Access.", width=130, type=["numericColumn", "numberColumnFilter"])

if 'Smiles' in df_display.columns:
    gb.configure_column("Smiles", hide=True)

grid_options = gb.build()

# Display grid
grid_response = AgGrid(
    df_display,
    gridOptions=grid_options,
    update_mode=GridUpdateMode.SELECTION_CHANGED,
    height=400,
    fit_columns_on_grid_load=False,
    allow_unsafe_jscode=True,
    theme='streamlit'
)

# --- STRUCTURE PREVIEW ---
selected_rows = grid_response['selected_rows']

if selected_rows is not None and len(selected_rows) > 0:
    st.subheader("Structure Preview")
    
    # Handle AgGrid return format (sometimes list, sometimes DataFrame)
    selected_row = selected_rows.iloc[0] if hasattr(selected_rows, 'iloc') else selected_rows[0]
    ligand_id = selected_row['Ligand']
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.markdown(f"**Ligand:** {ligand_id}")
        
        # --- SINGLE API CALL ---
        if 'Smiles' in selected_row and selected_row['Smiles']:
            # Fetch all data (Name, Link, IUPAC, Desc) in one go
            pc_data = fetch_pubchem_data(selected_row['Smiles'])
            
            if pc_data and pc_data['CommonName']:
                st.markdown(f"### {pc_data['CommonName']}")
                if pc_data['Link']:
                    st.markdown(f"🔗 [View on PubChem]({pc_data['Link']})")
                
                if pc_data['IUPAC']:
                    with st.expander("IUPAC Name"):
                        st.text(pc_data['IUPAC'])
                
                if pc_data['Description']:
                    with st.expander("Description"):
                        st.caption(pc_data['Description'])
            else:
                st.caption("No additional data found in PubChem.")
                
            # Show SMILES in expander to save space
            with st.expander("Show SMILES"):
                st.code(selected_row['Smiles'])

        # --- METRICS ---
        st.markdown(f"**Target Affinity:** {selected_row['Target (kcal/mol)']} kcal/mol")
        
        if 'Specificity' in selected_row and selected_row['Specificity'] is not None:
            st.markdown(f"**Specificity:** {selected_row['Specificity']}")
        else:
            st.markdown("**Specificity:** N/A")

        # --- LAZY POCKET RESIDUE CALCULATION ---
        st.markdown(f"**Interacting Residues:**")
        receptor_path = "data/interim/target_prep.pdbqt"
        ligand_path = f"data/results/target/{RUN_NAME}/poses/{ligand_id}_docked.pdbqt"
        
        if os.path.exists(receptor_path) and os.path.exists(ligand_path):
            with st.spinner("Analyzing interactions..."):
                residues = get_pocket_residues_from_pose(receptor_path, ligand_path, cutoff=3.5)
                st.info(residues)
        else:
            st.warning("Structure files not found.")
        
        # --- 3D VIEWER BUTTON ---
        st.markdown("")
        if st.button("Open in 3D Viewer", type="primary", key="open_3d_viewer"):
            st.session_state.py3dmol_preselected_ligand = ligand_id
            st.switch_page("pages/6_py3Dmol Visualization.py") # Ensure filename matches exactly
    
    with col2:
        if 'Smiles' in selected_row and selected_row['Smiles']:
            try:
                mol = Chem.MolFromSmiles(selected_row['Smiles'])
                if mol:
                    img = Draw.MolToImage(mol, size=(400, 300))
                    st.image(img, caption=f"Structure of {ligand_id}")
                else:
                    st.error("Invalid SMILES")
            except Exception:
                st.error("Could not render structure.")
        else:
            st.info("No SMILES data available.")


# --- DOWNLOAD BUTTON ---
csv_data = df_display.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Download Summary Report (CSV)",
    data=csv_data,
    file_name="summary_report.csv",
    mime="text/csv"
)
