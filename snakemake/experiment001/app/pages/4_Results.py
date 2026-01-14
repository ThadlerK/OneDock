import streamlit as st
import pandas as pd
import os
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
from rdkit import Chem
from rdkit.Chem import Draw

st.title("Results & Ranking")

result_file = "data/results/docking_report.csv"

if os.path.exists(result_file):
    df = pd.read_csv(result_file).sort_values(by="Affinity_kcal_mol")
    
    st.subheader("Docking Scores & Visualization")
    
    # Create two columns: Left for Table, Right for Structure
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Configure the interactive grid
        gb = GridOptionsBuilder.from_dataframe(df)
        gb.configure_selection('single', use_checkbox=False) # Allow single row selection
        gb.configure_grid_options(domLayout='normal')
        gridOptions = gb.build()

        # Display the grid and wait for user interaction
        grid_response = AgGrid(
            df,
            gridOptions=gridOptions,
            update_mode=GridUpdateMode.SELECTION_CHANGED,
            fit_columns_on_grid_load=True,
            height=400,
            theme='streamlit' # 'streamlit', 'alpine', 'balham'
        )

    with col2:
        st.write("### Structure Preview")
        
        # Check if a row is selected
        selected = grid_response['selected_rows']
        
        if selected is not None and len(selected) > 0:
            # Get SMILES from the selected row (AgGrid returns a list of dicts)
            # Note: AgGrid might wrap the row in a specific object depending on version
            # usually it is a dict like {'Receptor':..., 'Smiles':...}
            selected_row = selected.iloc[0] if isinstance(selected, pd.DataFrame) else pd.DataFrame(selected).iloc[0]
            
            smiles = selected_row.get("Smiles", "")
            
            if smiles and smiles != "unknown":
                st.code(smiles) # Show the string
                
                # Render the image using RDKit
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        img = Draw.MolToImage(mol, size=(300, 300))
                        st.image(img, caption=selected_row.get("Ligand", "Structure"))
                    else:
                        st.error("Invalid SMILES string")
                except Exception as e:
                    st.error(f"Could not render structure: {e}")
            else:
                st.info("No valid SMILES found for this entry.")
        else:
            st.info("👈 Click on a row to visualize the molecule.")
    
    st.subheader("Best Pose Visualization")
    st.write("Select a ligand to visualize (Placeholder logic)")
else:
    st.warning("No results found. Please run the docking pipeline first.")