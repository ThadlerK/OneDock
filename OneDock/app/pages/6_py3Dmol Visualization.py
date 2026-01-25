import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))
from utils import create_py3dmol_visualization, convert_pdbqt_to_pdb, analyze_protein_ligand_interactions

st.set_page_config(page_title = "py3Dmol Visualization")
st.title("3D Molecular Visualization with py3Dmol")
st.markdown("### Interactive Structure Viewer with py3Dmol")

# Display information about py3Dmol
st.info("""
**py3Dmol** enables an interactive 3D visualization of protein-ligand complexes.

Reference: Rego, N. and Koes, D. 3Dmol.js: molecular visualization with WebGL. 
Bioinformatics 31, 1322-1324 (2015)
""")

# Check if docking results exist
result_file = "data/results/docking_report_target.csv"
poses_dir = "data/results/target/poses"
receptor_file = "data/interim/target_prep.pdbqt"

if not os.path.exists(result_file):
    st.warning("No docking results found. Please run the docking pipeline first.")
    if st.button("Go to Docking"):
        st.switch_page("pages/2_Docking.py")
    st.stop()

# Load docking results
df = pd.read_csv(result_file)

# Load reference results if available for specificity calculation
ref_file = "data/results/docking_report_reference.csv"
if os.path.exists(ref_file):
    df_ref = pd.read_csv(ref_file)
    
    df = df.merge(
        df_ref[['Ligand', 'Affinity_kcal_mol']],
        on='Ligand',
        how='left',
        suffixes=('', '_ref')
    )
    df['Specificity'] = df['Affinity_kcal_mol'] - df['Affinity_kcal_mol_ref']
    has_specificity = True
else:
    has_specificity = False

df_sorted = df.sort_values(by="Affinity_kcal_mol")

# Load PoseBusters results if available
pb_results_file = "data/results/posebusters/posebusters_validation.csv"
has_pb_results = os.path.exists(pb_results_file)

if has_pb_results:
    pb_df = pd.read_csv(pb_results_file)
    df_sorted = df_sorted.merge(
        pb_df[['ligand_id', 'quality_score']], 
        left_on='Ligand', 
        right_on='ligand_id', 
        how='left'
    )
    df_sorted['quality_score'] = df_sorted['quality_score'].fillna(0)

st.subheader("Select Ligand for Visualization")

# --- INITIALIZE SESSION STATE ---
if 'affinity_cutoff' not in st.session_state:
    st.session_state.affinity_cutoff = -3.0
if 'specificity_cutoff' not in st.session_state:
    st.session_state.specificity_cutoff = 0.0
if 'py3dmol_selected_ligand' not in st.session_state:
    st.session_state.py3dmol_selected_ligand = None
if 'py3dmol_visualization' not in st.session_state:
    st.session_state.py3dmol_visualization = None
if 'py3dmol_just_created' not in st.session_state:
    st.session_state.py3dmol_just_created = False

# --- CUSTOM CUTOFF FILTERS ---
col1, col2 = st.columns(2)

with col1:
    affinity_cutoff = st.number_input(
        "Target Affinity ≤ (kcal/mol):",
        value=st.session_state.affinity_cutoff,
        step=0.5,
        key='py3d_affinity'
    )
    st.session_state.affinity_cutoff = affinity_cutoff

with col2:
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Specificity ≤:",
            value=st.session_state.specificity_cutoff,
            step=0.5,
            key='py3d_specificity'
        )
        st.session_state.specificity_cutoff = specificity_cutoff
    else:
        st.info("No reference data - specificity filter not available")
        specificity_cutoff = None

# Apply filters
if has_specificity and specificity_cutoff is not None:
    df_display = df_sorted[
        (df_sorted['Affinity_kcal_mol'] <= affinity_cutoff) &
        (df_sorted['Specificity'] <= specificity_cutoff)
    ]
else:
    df_display = df_sorted[
        df_sorted['Affinity_kcal_mol'] <= affinity_cutoff
    ]

# Display table with selection
st.write(f"Showing {len(df_display)} ligands:")

# Prepare display dataframe
display_cols = ['Ligand', 'Affinity_kcal_mol']
if has_pb_results:
    display_cols.append('quality_score')
if 'Smiles' in df_display.columns:
    display_cols.append('Smiles')

# Create format function for selectbox
def format_ligand_label(x):
    affinity = df_display[df_display['Ligand']==x]['Affinity_kcal_mol'].values[0]
    label = f"{x} (Affinity: {affinity:.2f}"
    
    if has_specificity and 'Specificity' in df_display.columns:
        specificity = df_display[df_display['Ligand']==x]['Specificity'].values[0]
        label += f", Specificity: {specificity:.2f}"
    
    label += " kcal/mol)"
    return label

# Check if a ligand was preselected from Summary page
if 'py3dmol_preselected_ligand' in st.session_state and st.session_state.py3dmol_preselected_ligand:
    preselected = st.session_state.py3dmol_preselected_ligand
    # Check if the preselected ligand is in the filtered list
    if preselected in df_display['Ligand'].tolist():
        default_index = df_display['Ligand'].tolist().index(preselected)
    else:
        default_index = 0
    # Clear the preselection
    st.session_state.py3dmol_preselected_ligand = None
else:
    default_index = 0

selected_ligand = st.selectbox(
    "Choose a ligand to visualize:",
    options=df_display['Ligand'].tolist(),
    format_func=format_ligand_label,
    index=default_index
)

# Store selected ligand in session state
st.session_state.py3dmol_selected_ligand = selected_ligand

# Reset the flag at page load (before checking for saved visualization)
if 'py3dmol_just_created' not in st.session_state:
    st.session_state.py3dmol_just_created = False
    
# Only reset to False if we're loading the page fresh (not after creating viz)
if st.session_state.py3dmol_just_created:
    # Mark that we've shown the newly created viz once
    st.session_state.py3dmol_just_created = False

# --- DISPLAY PREVIOUSLY SAVED VISUALIZATION ---
# Show saved visualization if it exists
if st.session_state.py3dmol_visualization is not None:
    viz_data = st.session_state.py3dmol_visualization
    
    st.subheader(f"3D Visualization: {viz_data['ligand_id']}")
    
    # Display the saved visualization
    components.html(viz_data['html_content'], height=650, scrolling=False)
    
    # Show ligand information
    ligand_info = viz_data['ligand_info']
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Binding Affinity", f"{ligand_info['Affinity_kcal_mol']:.2f} kcal/mol")
    
    with col2:
        if has_pb_results and 'quality_score' in ligand_info:
            quality = ligand_info['quality_score'] * 100
            st.metric("Quality Score", f"{quality:.1f}%")
        else:
            st.metric("Quality Score", "N/A")
    
    with col3:
        st.metric("Ligand ID", viz_data['ligand_id'])
    
    # Display settings used
    with st.expander("Visualization Settings Used"):
        settings = viz_data['settings']
        st.write(f"**Protein Style:** {settings['protein_style']}")
        st.write(f"**Ligand Style:** {settings['ligand_style']}")
        st.write(f"**Show Surface:** {settings['show_surface']}")
        st.write(f"**Show Interactions:** {settings['show_interactions']}")
        st.write(f"**Show Docking Box:** {settings['show_docking_box']}")
        st.write(f"**Show Pocket Residues:** {settings['show_pocket_residues']}")
    
    # Download button for saved visualization
    st.download_button(
        label="Download Visualization (HTML)",
        data=viz_data['html_content'],
        file_name=f"{viz_data['ligand_id']}_visualization.html",
        mime="text/html"
    )
    
    if st.button("Create New Visualization"):
        st.session_state.py3dmol_visualization = None
        st.rerun()
    
    st.stop()

# Visualization settings
st.subheader("Create New Visualization")

col1, col2, col3, col4 = st.columns(4)

with col1:
    protein_style = st.selectbox(
        "Protein Style:",
        ["cartoon", "line", "cross", "stick", "sphere", "ribbon", "trace"]
    )
    
with col2:
    ligand_style = st.selectbox(
        "Ligand Style:",
        ["stick", "sphere", "line", "cross", "ball_and_stick"]
    )

# Display options in 2x2 grid
col1, col2 = st.columns(2)
with col1:
    show_surface = st.checkbox("Show Surface", value=False)
    show_interactions = st.checkbox("Show Interacting Residues", value=True, 
                                     help="Highlights residues within 3.5Å of the ligand in yellow-green")
with col2:
    show_docking_box = st.checkbox("Show Docking Box", value=False)
    show_pocket_residues = st.checkbox("Show Binding Pocket", value=True,
                                        help="Shows defined pocket residues in cyan")

interaction_distance = 3.5

# Generate visualization
if st.button("Generate 3D Visualization", type="primary"):
    
    # Prepare file paths
    receptor_pdb = "data/interim/receptor_ready.pdb"
    ligand_pdbqt = os.path.join(poses_dir, f"{selected_ligand}_docked.pdbqt")
    ligand_pdb = f"data/results/py3dmol/{selected_ligand}_docked.pdb"
    
    # Ensure output directory exists
    os.makedirs("data/results/py3dmol", exist_ok=True)
    
    # Convert files if needed
    with st.spinner("Preparing files for visualization..."):
        # Convert receptor
        if not os.path.exists(receptor_pdb):
            if os.path.exists(receptor_file):
                convert_pdbqt_to_pdb(receptor_file, receptor_pdb)
            else:
                st.error("Receptor file not found")
                st.stop()
        
        # Convert ligand
        if os.path.exists(ligand_pdbqt):
            convert_pdbqt_to_pdb(ligand_pdbqt, ligand_pdb)
        else:
            st.error(f"Ligand file not found: {ligand_pdbqt}")
            st.stop()
    
    # Create visualization
    with st.spinner("Generating 3D visualization..."):
        # Load pocket residues from config
        from utils import load_config, calculate_box_from_residues
        config = load_config()
        pocket_residues_str = config.get('pocket_residues', '')
        pocket_residues = [int(r.strip()) for r in pocket_residues_str.split(',') if r.strip().isdigit()] if pocket_residues_str else None
        
        # Get docking box parameters from config if available
        box_center = None
        box_size = None
        if show_docking_box:
            center_x = config.get('center_x')
            center_y = config.get('center_y')
            center_z = config.get('center_z')
            size_x = config.get('size_x', 20)
            size_y = config.get('size_y', 20)
            size_z = config.get('size_z', 20)
            
            if all(v is not None for v in [center_x, center_y, center_z]):
                box_center = (float(center_x), float(center_y), float(center_z))
                box_size = (float(size_x), float(size_y), float(size_z))
            elif pocket_residues:
                # Calculate box from pocket residues if not in config
                original_receptor = "data/inputs/target.pdb"
                if os.path.exists(original_receptor):
                    box_center, box_size = calculate_box_from_residues(original_receptor, pocket_residues, padding=5.0)
        
        # Use original uploaded receptor file
        original_receptor = "data/inputs/target.pdb"
        
        # Only pass pocket_residues if the checkbox is checked
        pocket_to_display = pocket_residues if show_pocket_residues else None
        
        html_content = create_py3dmol_visualization(
            receptor_pdb=receptor_pdb,
            ligand_pdb=ligand_pdb,
            pocket_residues=pocket_to_display,
            original_receptor_pdb=original_receptor,
            protein_style=protein_style,
            ligand_style=ligand_style,
            show_surface=show_surface,
            show_interactions=show_interactions,
            show_docking_box=show_docking_box,
            box_center=box_center,
            box_size=box_size,
            width=900,
            height=700
        )
    
    if html_content:
        # Store visualization in session state and mark as just created
        st.session_state.py3dmol_visualization = {
            'html_content': html_content,
            'ligand_id': selected_ligand,
            'ligand_info': df_display[df_display['Ligand'] == selected_ligand].iloc[0].to_dict(),
            'settings': {
                'protein_style': protein_style,
                'ligand_style': ligand_style,
                'show_surface': show_surface,
                'show_interactions': show_interactions,
                'show_docking_box': show_docking_box,
                'show_pocket_residues': show_pocket_residues
            }
        }
        st.session_state.py3dmol_just_created = True
        
        st.subheader(f"3D Visualization: {selected_ligand}")
        
        # Display the interactive 3D viewer
        components.html(html_content, height=650, scrolling=False)
        
        # Show ligand information
        ligand_info = df_display[df_display['Ligand'] == selected_ligand].iloc[0]
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Binding Affinity", f"{ligand_info['Affinity_kcal_mol']:.2f} kcal/mol")
        
        with col2:
            if has_pb_results and 'quality_score' in ligand_info:
                quality = ligand_info['quality_score'] * 100
                st.metric("Quality Score", f"{quality:.1f}%")
            else:
                st.metric("Quality Score", "N/A")
        
        with col3:
            st.metric("Ligand ID", selected_ligand)
        
        # Interaction analysis
        if show_interactions:
            st.subheader("Interaction Analysis")
            
            with st.spinner("Analyzing protein-ligand interactions..."):
                interactions = analyze_protein_ligand_interactions(
                    receptor_pdb=receptor_pdb,
                    ligand_pdb=ligand_pdb,
                    distance_cutoff=interaction_distance
                )
            
            if interactions:
                # Display interaction summary
                interaction_types = {}
                for interaction in interactions:
                    itype = interaction.get('type', 'unknown')
                    interaction_types[itype] = interaction_types.get(itype, 0) + 1
                
                st.write("**Interaction Summary:**")
                cols = st.columns(len(interaction_types))
                for idx, (itype, count) in enumerate(interaction_types.items()):
                    with cols[idx]:
                        st.metric(itype.replace('_', ' ').title(), count)
                
                # Detailed interactions table
                interactions_df = pd.DataFrame(interactions)
                st.dataframe(interactions_df, use_container_width=True)
                
                # Download interactions
                csv = interactions_df.to_csv(index=False)
                st.download_button(
                    label="Download Interactions (CSV)",
                    data=csv,
                    file_name=f"{selected_ligand}_interactions.csv",
                    mime="text/csv"
                )
            else:
                st.info("No significant interactions detected")
        
        # Download visualization
        st.download_button(
            label="Download HTML Visualization",
            data=html_content,
            file_name=f"{selected_ligand}_visualization.html",
            mime="text/html"
        )
        
    else:
        st.error("Failed to generate visualization")
