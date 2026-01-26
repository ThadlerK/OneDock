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
def fetch_pubchem_info(smiles):
    """
    Queries PubChem PUG REST API to get the common name and URL for a SMILES string.
    """
    if not smiles or len(str(smiles)) < 2: 
        return None, None

    try:
        # Encode SMILES safely for URL
        encoded_smiles = urllib.parse.quote(smiles)
        
        # 1. Get the CID (Compound ID) from SMILES
        # /compound/smiles/{smiles}/cids/JSON
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        url_cid = f"{base_url}/compound/smiles/{encoded_smiles}/cids/JSON"
        
        resp = requests.get(url_cid, timeout=3)
        if resp.status_code == 200:
            cid = resp.json()['IdentifierList']['CID'][0]
            
            # 2. Get the Title (Common Name) using the CID
            # /compound/cid/{cid}/property/Title/JSON
            url_name = f"{base_url}/compound/cid/{cid}/property/Title/JSON"
            resp_name = requests.get(url_name, timeout=3)
            
            name = "Unknown"
            if resp_name.status_code == 200:
                props = resp_name.json()['PropertyTable']['Properties'][0]
                name = props.get('Title', f"PubChem CID: {cid}")
            
            # 3. Construct Link
            link = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            
            return name, link
            
    except Exception as e:
        return None, None
        
    return None, None

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
TARGET_FILE = f"data/results/docking_report_target_{lib_name}.csv"
REF_FILE = f"data/results/docking_report_reference_{lib_name}.csv"

@st.cache_data
def load_summary_data(target_path, ref_path):
    if not os.path.exists(target_path):
        return None
    
    # Load target docking results
    df = pd.read_csv(target_path)
    
    # Rename columns
    df = df.rename(columns={
        'Affinity_kcal_mol': 'Affinity_Target',
        'Ligand': 'Ligand'
    })
    
    # Round affinity to 2 decimal places
    df['Affinity_Target'] = df['Affinity_Target'].round(2)
    
    # Load reference if available and calculate specificity
    if os.path.exists(ref_path):
        df_ref = pd.read_csv(ref_path)
        df_ref = df_ref.rename(columns={
            'Affinity_kcal_mol': 'Affinity_Ref',
            'Ligand': 'Ligand'
        })
        
        df = df.merge(df_ref[['Ligand', 'Affinity_Ref']], on='Ligand', how='left')
        df['Specificity'] = (df['Affinity_Target'] - df['Affinity_Ref']).round(2)
    else:
        df['Specificity'] = None
    
    # Calculate pocket residues for each ligand
    receptor_pdbqt = "data/interim/target_prep.pdbqt"
    pocket_residues_list = []
    
    for idx, row in df.iterrows():
        ligand_id = row['Ligand']
        ligand_pdbqt = f"data/results/target/poses/{ligand_id}_docked.pdbqt"
        
        if os.path.exists(receptor_pdbqt) and os.path.exists(ligand_pdbqt):
            residues = get_pocket_residues_from_pose(receptor_pdbqt, ligand_pdbqt, cutoff=3.5)
            pocket_residues_list.append(residues)
        else:
            pocket_residues_list.append("N/A")
    
    df['Pocket_Residues'] = pocket_residues_list
    
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

    if results_files:
        mmpbsa_data = []
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
                                "Delta G (MMPBSA)": val
                            })
                            break
            except:
                pass
        
    if mmpbsa_data and df is not None:
        df_mmpbsa = pd.DataFrame(mmpbsa_data)
        # B. Handle multiple ranks (Optional but recommended)
        df_mmpbsa_best = df_mmpbsa.sort_values("Delta G (MMPBSA)").groupby("Ligand").first().reset_index()

        # C. Merge into main dataframe
        df = df.merge(df_mmpbsa_best[['Ligand', 'Delta G (MMPBSA)']], on='Ligand', how='left')
                    
    
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
    st.session_state.affinity_cutoff = -3.0
if 'specificity_cutoff' not in st.session_state:
    st.session_state.specificity_cutoff = 0.0

# --- FILTER CONTROLS ---
has_specificity = df['Specificity'].notna().any()

col1, col2 = st.columns(2)

with col1:
    affinity_cutoff = st.number_input(
        "Target Affinity ≤ (kcal/mol):",
        value=st.session_state.affinity_cutoff,
        step=0.5,
        key='summary_affinity'
    )
    st.session_state.affinity_cutoff = affinity_cutoff

with col2:
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Specificity ≤:",
            value=st.session_state.specificity_cutoff,
            step=0.5,
            key='summary_specificity'
        )
        st.session_state.specificity_cutoff = specificity_cutoff
    else:
        st.info("No reference data - specificity filter not available")
        specificity_cutoff = None

# --- APPLY FILTERS ---
if has_specificity and specificity_cutoff is not None:
    df = df[
        (df['Affinity_Target'] <= affinity_cutoff) &
        (df['Specificity'] <= specificity_cutoff)
    ]
else:
    df = df[df['Affinity_Target'] <= affinity_cutoff]

# --- DISPLAY SUMMARY TABLE ---
st.subheader("Summary Table")

# Prepare display dataframe
display_columns = ['Ligand', 'Affinity_Target', 'Specificity', 'Pocket_Residues']
if 'PoseBusters_Score' in df.columns:
    display_columns.append('PoseBusters_Score')
if 'Lipinski_RO5' in df.columns:
    display_columns.append('Lipinski_RO5')
if 'Smiles' in df.columns:
    display_columns.append('Smiles')
if 'Delta G (MMPBSA)' in df.columns:
    display_columns.append('Delta G (MMPBSA)')

df_display = df[display_columns].copy()

# Rename for better display
rename_dict = {
    'Affinity_Target': 'Target (kcal/mol)',
    'Pocket_Residues': 'Pocket Residues'
}
if 'PoseBusters_Score' in df_display.columns:
    rename_dict['PoseBusters_Score'] = 'PoseBusters (%)'
if 'Lipinski_RO5' in df_display.columns:
    rename_dict['Lipinski_RO5'] = 'Lipinski RO5'

df_display = df_display.rename(columns=rename_dict)

# Configure AgGrid
gb = GridOptionsBuilder.from_dataframe(df_display)
gb.configure_selection(selection_mode='single', use_checkbox=False)
gb.configure_column("Ligand", width=120)
gb.configure_column("Target (kcal/mol)", width=150, type=["numericColumn", "numberColumnFilter"])
gb.configure_column("Specificity", width=120, type=["numericColumn", "numberColumnFilter"])
gb.configure_column("Pocket Residues", width=400, wrapText=True, autoHeight=True)
if 'PoseBusters (%)' in df_display.columns:
    gb.configure_column("PoseBusters (%)", width=130, type=["numericColumn", "numberColumnFilter"])
if 'Lipinski RO5' in df_display.columns:
    gb.configure_column("Lipinski RO5", width=150)

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
    
    selected_row = selected_rows.iloc[0] if hasattr(selected_rows, 'iloc') else selected_rows[0]
    ligand_id = selected_row['Ligand']
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.markdown(f"**Ligand:** {ligand_id}")

        if 'Smiles' in selected_row and selected_row['Smiles']:
            real_name, pubchem_link = fetch_pubchem_info(selected_row['Smiles'])
            
            if real_name:
                st.markdown(f"### {real_name}") # Display Name Big
                st.markdown(f"[View on PubChem]({pubchem_link})")
            else:
                st.caption("Name not found in PubChem")
        
        st.markdown(f"**Target Affinity:** {selected_row['Target (kcal/mol)']} kcal/mol")
        
        if 'Specificity' in selected_row and selected_row['Specificity'] is not None:
            st.markdown(f"**Specificity:** {selected_row['Specificity']}")
        else:
            st.markdown("**Specificity:** N/A (no reference)")
        
        st.markdown(f"**Pocket Residues:**")
        st.text(selected_row['Pocket Residues'])
        
        if 'Smiles' in selected_row:
            st.markdown(f"**SMILES:** `{selected_row['Smiles']}`")
        
        # Add button to open in 3D Viewer
        st.markdown("")
        if st.button("Open in 3D Viewer", type="primary", key="open_3d_viewer"):
            # Store the selected ligand in session state for py3Dmol to use
            st.session_state.py3dmol_preselected_ligand = ligand_id
            st.switch_page("pages/6_py3Dmol Visualization.py")
    
    with col2:
        if 'Smiles' in selected_row and selected_row['Smiles']:
            try:
                mol = Chem.MolFromSmiles(selected_row['Smiles'])
                if mol:
                    img = Draw.MolToImage(mol, size=(400, 300))
                    st.image(img, caption=f"Structure of {ligand_id}")
                else:
                    st.error("Could not generate structure from SMILES")
            except Exception as e:
                st.error(f"Error generating structure: {e}")
        else:
            st.info("No SMILES data available for structure preview")

# --- DOWNLOAD BUTTON ---
csv_data = df_display.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Download Summary Report (CSV)",
    data=csv_data,
    file_name="summary_report.csv",
    mime="text/csv"
)
