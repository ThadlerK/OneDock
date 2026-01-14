#!/bin/bash

# --- CONFIGURATION ---
BASE_DIR="/workspace/vina_experiment004"

# 1. RECEPTOR SELECTION
# Set to ("all") to dock against every PDB in the structures directory.
# OR list specific filenames (without .pdb extension) to limit the run.
# Example: SELECTED_RECEPTORS=("SLC6A20_8I91_occluded" "SLC6A20_inward")
SELECTED_RECEPTORS=("SLC6A20_8WM3_complete_inwardopen")

# 2. LIGAND LIBRARY SELECTION (NEW)
# Set to ("all") to dock every .sml file found in data/library1_smiles
# OR list specific filenames (without .sml extension)
# Example: SELECTED_LIBRARIES=("Fragment_Library" "ChEMBL_Drugs")
SELECTED_LIBRARIES=("all")

# 3. DOCKING TARGET CONFIGURATION
TARGET_SITE="orthosteric"   # Options: "orthosteric", "allosteric", "all"
TARGET_CONF="auto"          # Options: "occluded", "inward", "outward", "auto" (detects from filename)

# 4. PATHS
BINDING_DATA_FILE="$BASE_DIR/data/reference/pockets.csv"
STRUCTURES_DIR="$BASE_DIR/data/structures"
LIGAND_SML_DIR="$BASE_DIR/data/library1_smiles"

# Intermediate/Output Directories
LIGAND_PDBQT_DIR="$BASE_DIR/data/library1_pdbqt"
OUTPUT_BASE_DIR="$BASE_DIR/outputs/${TARGET_SITE}_${TARGET_CONF}"

# 5. VINA SETTINGS
GRID_SIZE=20          
EXHAUSTIVENESS=8      

# --- SETUP ---
mkdir -p "$OUTPUT_BASE_DIR"
PYTHON_SCRIPT="$BASE_DIR/scripts/get_pocket_center.py"

echo "========================================"
echo "STARTING DOCKING PIPELINE"
echo "Target Site: $TARGET_SITE"
echo "Target Conf: $TARGET_CONF"
echo "Receptor Filter: ${SELECTED_RECEPTORS[*]}"
echo "Library Filter:  ${SELECTED_LIBRARIES[*]}"
echo "========================================"

# --- OUTER LOOP: RECEPTORS ---
for RECEPTOR_PDB in "$STRUCTURES_DIR"/*.pdb; do
    
    REC_FILENAME=$(basename -- "$RECEPTOR_PDB")
    REC_NAME="${REC_FILENAME%.*}"
    RECEPTOR_PDBQT="$STRUCTURES_DIR/${REC_NAME}.pdbqt"

    # --- RECEPTOR SELECTION LOGIC ---
    PROCESS_THIS_REC=false
    if [ "${SELECTED_RECEPTORS[0]}" == "all" ]; then
        PROCESS_THIS_REC=true
    else
        for TARGET in "${SELECTED_RECEPTORS[@]}"; do
            if [[ "$REC_NAME" == "$TARGET" ]]; then
                PROCESS_THIS_REC=true; break
            fi
        done
    fi

    if [ "$PROCESS_THIS_REC" = false ]; then continue; fi
    # -------------------------------

    # 1. CALCULATE CENTER
    COORDS=$(python3 "$PYTHON_SCRIPT" "$RECEPTOR_PDB" "$BINDING_DATA_FILE" "$REC_NAME" --site "$TARGET_SITE" --conf "$TARGET_CONF")

    if [ -z "$COORDS" ]; then
        echo "  > Skipping $REC_NAME: No binding site found."
        continue
    fi

    read CENTER_X CENTER_Y CENTER_Z <<< "$COORDS"
    echo "  > Processing $REC_NAME (Pocket: $CENTER_X $CENTER_Y $CENTER_Z)"

    # --- CREATE RECEPTOR SPECIFIC FOLDER ---
    REC_OUTPUT_DIR="$OUTPUT_BASE_DIR/$REC_NAME"
    OUTPUT_POSE_DIR="$REC_OUTPUT_DIR/pdbqt"
    LOG_FILE_DIR="$REC_OUTPUT_DIR/logs"
    mkdir -p "$REC_OUTPUT_DIR"
    mkdir -p "$OUTPUT_POSE_DIR"
    mkdir -p "$LOG_FILE_DIR"

    # --- CREATE RECEPTOR SPECIFIC SUMMARY FILE ---
    SUMMARY_FILE="$REC_OUTPUT_DIR/${REC_NAME}_summary.csv"
    if [ ! -f "$SUMMARY_FILE" ]; then
        echo "Receptor,Ligand_Library,Ligand_ID,Affinity_kcal_mol,SMILES" > "$SUMMARY_FILE"
    fi

    # Convert Receptor to PDBQT if needed
    if [ ! -f "$RECEPTOR_PDBQT" ]; then
        obabel -ipdb "$RECEPTOR_PDB" -opdbqt -O "$RECEPTOR_PDBQT" -xr -h --partialcharge gasteiger 2> /dev/null
    fi

    # --- MIDDLE LOOP: SML FILES ---
    for LIGAND_SML in "$LIGAND_SML_DIR"/*.sml; do
        
        LIG_FILENAME=$(basename -- "$LIGAND_SML")
        LIG_NAME="${LIG_FILENAME%.*}" 

        # --- LIBRARY SELECTION LOGIC ---
        PROCESS_THIS_LIB=false
        if [ "${SELECTED_LIBRARIES[0]}" == "all" ]; then
            PROCESS_THIS_LIB=true
        else
            for TARGET_LIB in "${SELECTED_LIBRARIES[@]}"; do
                if [[ "$LIG_NAME" == "$TARGET_LIB" ]]; then
                    PROCESS_THIS_LIB=true; break
                fi
            done
        fi

        if [ "$PROCESS_THIS_LIB" = false ]; then continue; fi
        # -------------------------------

        echo "    ------------------------------------"
        echo "    Processing Library: $LIG_NAME"

        # 2. PYTHON: SPLIT & CONVERT
        python3 -c "
import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

sml_file = sys.argv[1]
out_dir = sys.argv[2]
base_name = sys.argv[3]

if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

try:
    with open(sml_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]

    print(f'      > Found {len(lines)} SMILES in file.')

    for idx, line in enumerate(lines):
        smiles = line.split()[0]
        
        # Use 5-digit padding as requested
        out_filename = f'{base_name}_{idx:05d}.pdbqt'
        smi_filename = f'{base_name}_{idx:05d}.smi'
        
        out_path = os.path.join(out_dir, out_filename)
        smi_path = os.path.join(out_dir, smi_filename)

        if os.path.exists(out_path): continue

        # Save SMILES text
        with open(smi_path, 'w') as f_smi: f_smi.write(smiles)

        mol = Chem.MolFromSmiles(smiles)
        if not mol: continue

        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            pdbqt_data = PDBQTWriterLegacy.write_string(mol_setups[0])
            pdbqt_string = pdbqt_data[0] if isinstance(pdbqt_data, tuple) else pdbqt_data
            with open(out_path, 'w') as out_f: out_f.write(pdbqt_string)
        except: pass
except: sys.exit(1)
" "$LIGAND_SML" "$LIGAND_PDBQT_DIR" "$LIG_NAME"

        # --- INNER LOOP: DOCKING ---
        shopt -s nullglob
        PDBQT_FILES=("$LIGAND_PDBQT_DIR/${LIG_NAME}_"*.pdbqt)
        shopt -u nullglob

        for SUB_LIGAND_PATH in "${PDBQT_FILES[@]}"; do
            
            SUB_FILENAME=$(basename -- "$SUB_LIGAND_PATH")
            SUB_LIG_ID="${SUB_FILENAME%.*}" 
            
            OUTPUT_POSE="$OUTPUT_POSE_DIR/${SUB_LIG_ID}_docked.pdbqt"
            LOG_FILE="$LOG_FILE_DIR/${SUB_LIG_ID}_docked.log"

            # Check specific summary file
            if [ -f "$SUMMARY_FILE" ] && grep -q "$REC_NAME,$LIG_NAME,$SUB_LIG_ID" "$SUMMARY_FILE"; then
                 continue
            fi

            echo "      > Docking: $SUB_LIG_ID"

            vina \
                --receptor "$RECEPTOR_PDBQT" \
                --ligand "$SUB_LIGAND_PATH" \
                --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
                --size_x $GRID_SIZE --size_y $GRID_SIZE --size_z $GRID_SIZE \
                --exhaustiveness $EXHAUSTIVENESS \
                --out "$OUTPUT_POSE" \
                > "$LOG_FILE" 2>&1

            if [ -f "$LOG_FILE" ]; then
                AFFINITY=$(grep "   1 " "$LOG_FILE" | head -n 1 | awk '{print $2}')
                
                SMI_FILE="$LIGAND_PDBQT_DIR/${SUB_LIG_ID}.smi"
                if [ -f "$SMI_FILE" ]; then CURRENT_SMILES=$(cat "$SMI_FILE"); else CURRENT_SMILES="NA"; fi
                
                # Write to Receptor-Specific Summary File
                echo "$REC_NAME,$LIG_NAME,$SUB_LIG_ID,$AFFINITY,$CURRENT_SMILES" >> "$SUMMARY_FILE"
            fi
        
        done
    done
done
echo "Batch Docking Complete."