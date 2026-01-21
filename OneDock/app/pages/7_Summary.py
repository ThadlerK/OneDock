import streamlit as st
import pandas as pd
import numpy as np
import os
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(layout="wide", page_title="Summary Report")

st.title("Summary Report")

# --- HELPER FUNCTIONS ---

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
TARGET_FILE = "data/results/docking_report_target.csv"
REF_FILE = "data/results/docking_report_reference.csv"

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
    
    # Sort by affinity (best first)
    df = df.sort_values('Affinity_Target', ascending=True)
    
    return df

df = load_summary_data(TARGET_FILE, REF_FILE)

if df is None:
    st.warning("No docking results found. Please run the docking pipeline first.")
    if st.button("Go to Docking"):
        st.switch_page("pages/2_Docking.py")
    st.stop()

# --- DISPLAY SUMMARY TABLE ---
st.subheader("Summary Table")

# Prepare display dataframe
display_columns = ['Ligand', 'Affinity_Target', 'Specificity', 'Pocket_Residues']
if 'Smiles' in df.columns:
    display_columns.append('Smiles')

df_display = df[display_columns].copy()

# Rename for better display
df_display = df_display.rename(columns={
    'Affinity_Target': 'Target (kcal/mol)',
    'Pocket_Residues': 'Pocket Residues'
})

# Configure AgGrid
gb = GridOptionsBuilder.from_dataframe(df_display)
gb.configure_selection(selection_mode='single', use_checkbox=False)
gb.configure_column("Ligand", width=120)
gb.configure_column("Target (kcal/mol)", width=150, type=["numericColumn", "numberColumnFilter"])
gb.configure_column("Specificity", width=120, type=["numericColumn", "numberColumnFilter"])
gb.configure_column("Pocket Residues", width=400, wrapText=True, autoHeight=True)

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
        st.markdown(f"**Target Affinity:** {selected_row['Target (kcal/mol)']} kcal/mol")
        
        if 'Specificity' in selected_row and selected_row['Specificity'] is not None:
            st.markdown(f"**Specificity:** {selected_row['Specificity']}")
        else:
            st.markdown("**Specificity:** N/A (no reference)")
        
        st.markdown(f"**Pocket Residues:**")
        st.text(selected_row['Pocket Residues'])
        
        if 'Smiles' in selected_row:
            st.markdown(f"**SMILES:** `{selected_row['Smiles']}`")
    
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
