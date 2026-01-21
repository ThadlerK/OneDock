import streamlit as st
import pandas as pd
import os
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent))
from utils import run_posebusters_validation, convert_pdbqt_to_pdb

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
result_file = "data/results/docking_report_target.csv"
poses_dir = "data/results/target/poses"

if not os.path.exists(result_file):
    st.warning("No docking results found. Please run the docking pipeline first.")
    if st.button("Go to Docking"):
        st.switch_page("pages/3_Docking.py")
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

st.subheader("PoseBusters Configuration")

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
        key='pb_affinity'
    )
    st.session_state.affinity_cutoff = affinity_cutoff

with col2:
    if has_specificity:
        specificity_cutoff = st.number_input(
            "Specificity ≤:",
            value=st.session_state.specificity_cutoff,
            step=0.5,
            key='pb_specificity'
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
    ]
else:
    selected_df = df_sorted[
        df_sorted['Affinity_kcal_mol'] <= affinity_cutoff
    ]

selected_ligands = selected_df['Ligand'].tolist()

# Display selection
display_cols = ['Ligand', 'Affinity_kcal_mol']
if has_specificity:
    display_cols.append('Specificity')

st.dataframe(selected_df[display_cols], height=250)

# Run PoseBusters
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
