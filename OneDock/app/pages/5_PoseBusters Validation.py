import streamlit as st
import pandas as pd
import os
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))
from utils import run_posebusters_validation, convert_pdbqt_to_pdb, load_config

st.set_page_config(page_title = "PoseBusters")
st.title("PoseBusters Validation")
st.markdown("### Quality Assessment of Docking Poses")

# Display information about PoseBusters
st.info("""
**PoseBusters** validates the physical and chemical plausibility of protein-ligand complexes.

Quality Score Interpretation:

:green[80-100%: High quality pose]

:orange[70-79%: Moderate quality]

:red[<70%: Low quality]

Reference: Buttenschoen, M., Morris, G.M. & Deane, C.M. 
PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. 
Chem. Sci. 15, 3130-3139 (2024)
""")

# Check if docking results exist
config = load_config()

lib_name = config.get("library_name", '')
result_file = f"data/results/docking_report_target_{lib_name}.csv"
poses_dir = "data/results/target/poses"

if not os.path.exists(result_file):
    st.warning("No docking results found. Please run the docking pipeline first.")
    if st.button("Go to Docking"):
        st.switch_page("pages/2_Docking.py")
    st.stop()

# Load docking results
df = pd.read_csv(result_file)

# Load reference results if available for specificity and rank gain calculation
ref_file = "data/results/docking_report_reference.csv"
if os.path.exists(ref_file):
    df_ref = pd.read_csv(ref_file)
    
    # Sort by affinity to assign ranks
    df_target_sorted = df.sort_values('Affinity_kcal_mol').reset_index(drop=True)
    df_target_sorted['Rank_Target'] = df_target_sorted.index + 1
    
    df_ref_sorted = df_ref.sort_values('Affinity_kcal_mol').reset_index(drop=True)
    df_ref_sorted['Rank_Ref'] = df_ref_sorted.index + 1
    
    # Merge ranks and reference affinity
    df = df.merge(df_target_sorted[['Ligand', 'Rank_Target']], on='Ligand', how='left')
    df = df.merge(df_ref_sorted[['Ligand', 'Rank_Ref']], on='Ligand', how='left')
    df = df.merge(
        df_ref[['Ligand', 'Affinity_kcal_mol']],
        on='Ligand',
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

st.subheader("PoseBusters Configuration")

# --- INITIALIZE SESSION STATE ---
if 'affinity_cutoff' not in st.session_state:
    st.session_state.affinity_cutoff = -3.0
if 'specificity_cutoff' not in st.session_state:
    st.session_state.specificity_cutoff = 0.0
if 'posebusters_results' not in st.session_state:
    st.session_state.posebusters_results = None
if 'pb_selected_ligands' not in st.session_state:
    st.session_state.pb_selected_ligands = None

# --- CUSTOM CUTOFF FILTERS ---
col1, col2, col3 = st.columns(3)

with col1:
    affinity_cutoff = st.number_input(
        "Max Target Affinity (kcal/mol):",
        value=st.session_state.affinity_cutoff,
        step=0.5,
        key='pb_affinity',
        help="Show only ligands with affinity ≤ this value."
    )
    st.session_state.affinity_cutoff = affinity_cutoff

with col2:
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Max Specificity:",
            value=st.session_state.specificity_cutoff,
            step=0.5,
            key='pb_specificity',
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
            key='pb_rank_gain',
            help="Show only ligands with rank gain ≥ this value."
        )
        st.session_state.rank_gain_cutoff = rank_gain_cutoff
    else:
        st.warning("⚠️ No reference data - rank gain filter not available")
        rank_gain_cutoff = None

# Apply filters
selected_df = df_sorted[df_sorted['Affinity_kcal_mol'] <= affinity_cutoff]

if has_specificity and specificity_cutoff is not None:
    selected_df = selected_df[selected_df['Specificity'] <= specificity_cutoff]

if has_rank_gain and rank_gain_cutoff is not None:
    selected_df = selected_df[selected_df['Rank_Gain'] >= rank_gain_cutoff]

selected_ligands = selected_df['Ligand'].tolist()

# Display selection
display_cols = ['Ligand', 'Affinity_kcal_mol']
if has_specificity:
    display_cols.append('Specificity')
if has_rank_gain:
    display_cols.append('Rank_Gain')

st.dataframe(selected_df[display_cols], height=250)

# --- DISPLAY PREVIOUSLY SAVED RESULTS ---
if st.session_state.posebusters_results is not None:
    st.markdown("---")
    st.subheader("Previously Calculated PoseBusters Results")
    st.info("These results were saved from a previous validation and persist across page navigation")
    
    results_df = st.session_state.posebusters_results
    
    # Display summary statistics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        avg_score = results_df['quality_score'].mean() * 100
        st.metric("Average Quality Score", f"{avg_score:.1f}%")
    
    with col2:
        passed = (results_df['quality_score'] >= 0.8).sum()
        st.metric("Poses Passed (>80%)", f"{passed}/{len(results_df)}")
    
    with col3:
        failed = (results_df['quality_score'] < 0.7).sum()
        st.metric("Poses Failed (<70%)", failed)
    
    with col4:
        best_score = results_df['quality_score'].max() * 100
        st.metric("Best Score", f"{best_score:.1f}%")
    
    # Display results table
    display_df = results_df.copy()
    display_df['quality_score'] = (display_df['quality_score'] * 100).round(1)
    st.dataframe(display_df, height=300)
    
    # Download button
    csv = results_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Stored PoseBusters Results",
        data=csv,
        file_name="stored_posebusters_results.csv",
        mime="text/csv"
    )
    
    if st.button("Clear Stored PoseBusters Results"):
        st.session_state.posebusters_results = None
        st.session_state.pb_selected_ligands = None
        st.rerun()

# Run PoseBusters
st.markdown("---")
if st.button("Run PoseBusters Validation", type="primary"):
    if not selected_ligands:
        st.error("Please select at least one ligand for validation")
        st.stop()
    
    progress_bar = st.progress(0)
    
    # Prepare output directory
    pb_output_dir = "data/results/posebusters"
    os.makedirs(pb_output_dir, exist_ok=True)
    
    # Convert PDBQT to PDB for selected ligands
    pdb_files = []
    receptor_pdb = "data/interim/receptor_ready.pdb"
    
    # Convert receptor if needed
    if not os.path.exists(receptor_pdb):
        receptor_pdbqt = "data/interim/receptor_ready.pdbqt"
        if os.path.exists(receptor_pdbqt):
            with st.spinner("Converting receptor to PDB..."):
                convert_pdbqt_to_pdb(receptor_pdbqt, receptor_pdb)
    
    # Convert ligands
    for idx, ligand_id in enumerate(selected_ligands):
        progress = (idx + 1) / len(selected_ligands)
        progress_bar.progress(progress)
        
        pdbqt_file = os.path.join(poses_dir, f"{ligand_id}_docked.pdbqt")
        pdb_file = os.path.join(pb_output_dir, f"{ligand_id}_docked.pdb")
        
        if os.path.exists(pdbqt_file):
            convert_pdbqt_to_pdb(pdbqt_file, pdb_file)
            pdb_files.append((ligand_id, pdb_file))
    
    progress_bar.empty()
    
    # Run PoseBusters validation
    with st.spinner("Analyzing poses..."):
        results_df = run_posebusters_validation(
            receptor_pdb=receptor_pdb,
            ligand_pdb_files=pdb_files,
            output_dir=pb_output_dir
        )
    
    if results_df is not None and not results_df.empty:
        # Store in session state
        st.session_state.posebusters_results = results_df
        st.session_state.pb_selected_ligands = selected_ligands
        
        # Display summary statistics
        st.subheader("Validation Summary")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            avg_score = results_df['quality_score'].mean() * 100
            st.metric("Average Quality Score", f"{avg_score:.1f}%")
        
        with col2:
            passed = (results_df['quality_score'] >= 0.8).sum()
            st.metric("Poses Passed (>80%)", f"{passed}/{len(results_df)}")
        
        with col3:
            failed = (results_df['quality_score'] < 0.7).sum()
            st.metric("Poses Failed (<70%)", failed)
        
        with col4:
            best_score = results_df['quality_score'].max() * 100
            st.metric("Best Score", f"{best_score:.1f}%")
        
        # Display detailed results
        st.subheader("Detailed Results")
        
        # Add color coding for quality scores and individual test results
        def color_quality_score(val):
            if val >= 80:
                return 'background-color: #90EE90'  # Light green
            elif val >= 70:
                return 'background-color: #FFD700'  # Gold
            else:
                return 'background-color: #FFB6C6'  # Light red
        
        def color_test_result(val):
            """No coloring for test results"""
            return ''  # No color
        
        # Format the dataframe and reorder columns
        display_df = results_df.copy()
        display_df['quality_score'] = display_df['quality_score'] * 100
        display_df = display_df.round(2)
        
        # Define column order matching Test Descriptions and rename with descriptions
        column_order = [
            'ligand_id',
            'quality_score',
            'mol_pred_loaded',
            'mol_cond_loaded',
            'sanitization',
            'inchi_convertible',
            'all_atoms_connected',
            'no_radicals',
            'bond_lengths',
            'bond_angles',
            'internal_steric_clash',
            'aromatic_ring_flatness',
            'non-aromatic_ring_non-flatness',
            'double_bond_flatness',
            'internal_energy',
            'minimum_distance_to_protein',
            'volume_overlap_with_protein',
            'minimum_distance_to_organic_cofactors',
            'volume_overlap_with_organic_cofactors',
            'minimum_distance_to_inorganic_cofactors',
            'volume_overlap_with_inorganic_cofactors',
            'minimum_distance_to_waters',
            'volume_overlap_with_waters',
            'passed_tests',
            'total_tests',
            'failed_tests_summary'
        ]
        
        # Filter to only include existing columns
        column_order = [col for col in column_order if col in display_df.columns]
        display_df = display_df[column_order]
        
        # Rename columns with descriptions matching Test Descriptions
        column_descriptions = {
            'ligand_id': 'Ligand ID',
            'quality_score': 'Quality Score (%)',
            'mol_pred_loaded': 'Mol Pred Loaded',
            'mol_cond_loaded': 'Mol Cond Loaded',
            'sanitization': 'Sanitization',
            'inchi_convertible': 'InChI Convertible',
            'all_atoms_connected': 'All Atoms Connected',
            'no_radicals': 'No Radicals',
            'bond_lengths': 'Bond Lengths Valid',
            'bond_angles': 'Bond Angles Valid',
            'internal_steric_clash': 'No Internal Clashes',
            'aromatic_ring_flatness': 'Aromatic Rings Flat',
            'non-aromatic_ring_non-flatness': 'Non-Aromatic Ring Non-Flatness',
            'double_bond_flatness': 'Double Bonds Flat',
            'internal_energy': 'Internal Energy OK',
            'minimum_distance_to_protein': 'Min Distance to Protein',
            'volume_overlap_with_protein': 'No Volume Overlap (Protein)',
            'minimum_distance_to_organic_cofactors': 'Min Distance to Cofactors',
            'volume_overlap_with_organic_cofactors': 'No Volume Overlap (Org. Cofactors)',
            'minimum_distance_to_inorganic_cofactors': 'Min Distance to Inorg. Cofactors',
            'volume_overlap_with_inorganic_cofactors': 'No Volume Overlap (Inorg. Cofactors)',
            'minimum_distance_to_waters': 'Min Distance to Waters',
            'volume_overlap_with_waters': 'No Volume Overlap (Waters)',
            'passed_tests': 'Passed Tests',
            'total_tests': 'Total Tests',
            'failed_tests_summary': 'Failed Tests Summary'
        }
        
        display_df = display_df.rename(columns=column_descriptions)
        
        # Apply styling: quality score gets color based on score, test results get pass/fail colors
        test_columns = [col for col in display_df.columns if col not in ['Ligand ID', 'Quality Score (%)', 'Passed Tests', 'Total Tests', 'Failed Tests Summary']]
        
        styled_df = display_df.style.applymap(
            color_quality_score,
            subset=['Quality Score (%)']
        ).applymap(
            color_test_result,
            subset=test_columns
        )
        
        st.dataframe(styled_df, use_container_width=True, height=400)
        
        # Show test descriptions
        with st.expander("Parameter Information"):
            st.markdown("""
            **All Atoms Connected:** Checks if all atoms in the ligand are part of a single connected molecule.
            
            **Bond Lengths Valid:** Validates that all bond lengths are within acceptable ranges.
            
            **Bond Angles Valid:** Checks that bond angles are chemically reasonable.
            
            **No Internal Clashes:** Ensures no steric clashes within the ligand itself.
            
            **Aromatic Rings Flat:** Verifies that aromatic ring systems are planar.
            
            **Double Bonds Flat:** Checks planarity of double bonds.
            
            **Internal Energy OK:** Validates that the internal energy of the ligand is reasonable.
            
            **Min Distance to Protein:** Ensures minimum distance between ligand and protein atoms.
            
            **No Volume Overlap (Protein):** Checks for volumetric overlaps with protein.
            
            **Min Distance to Cofactors:** Ensures minimum distance to organic cofactors.
            
            **No Volume Overlap (Org. Cofactors):** Checks for volumetric overlaps with organic cofactors.
            
            **Min Distance to Inorg. Cofactors:** Ensures minimum distance to inorganic cofactors.
            
            **No Volume Overlap (Inorg. Cofactors):** Checks for volumetric overlaps with inorganic cofactors.
            
            **Min Distance to Waters:** Ensures minimum distance to water molecules.
            
            **No Volume Overlap (Waters):** Checks for volumetric overlaps with waters.
            """)
        
        # Download results
        csv = results_df.to_csv(index=False)
        st.download_button(
            label="Download PoseBusters Results (CSV)",
            data=csv,
            file_name="posebusters_validation.csv",
            mime="text/csv"
        )
        
    else:
        st.error("PoseBusters validation failed. Please check the logs.")
