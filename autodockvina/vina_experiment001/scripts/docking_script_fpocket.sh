#!/bin/bash

# --- CONFIGURATION ---
BASE_DIR="/workspace/vina_experiment001"

STRUCTURES_DIR="$BASE_DIR/data/structures"
LIGAND_SML_DIR="$BASE_DIR/data/library1_smiles"
LIGAND_PDBQT_DIR="$BASE_DIR/data/library1_pdbqt"
OUTPUT_BASE_DIR="$BASE_DIR/outputs"

# Vina Grid Size for Pockets (Usually smaller than blind docking)
SIZE_X=18
SIZE_Y=18
SIZE_Z=18
EXHAUSTIVENESS=8

# --- PRE-FLIGHT CHECKS ---
echo "=== Phase 0: Preparing Environment ==="
mkdir -p "$LIGAND_PDBQT_DIR"
mkdir -p "$OUTPUT_BASE_DIR"

if ! command -v obabel &> /dev/null; then echo "Error: obabel not found"; exit 1; fi

# --- PHASE 1: PREPARE ALL LIGANDS (Run Once) ---
echo "=== Phase 1: Preparing Ligand Library ==="

for LIGAND_SML in "$LIGAND_SML_DIR"/*.sml; do
    LIG_FILENAME=$(basename -- "$LIGAND_SML")
    LIG_NAME="${LIG_FILENAME%.*}" # e.g., "Glutamine"
    
    # We pass the Input File, Output Directory, and Base Name to Python
    echo "  Processing library file: $LIG_NAME..."
    
    python3 -c "
import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

sml_file = sys.argv[1]
out_dir = sys.argv[2]
base_name = sys.argv[3]

try:
    with open(sml_file, 'r') as f:
        lines = f.readlines()

    count = 0
    for i, line in enumerate(lines):
        smiles = line.strip()
        if not smiles: continue # Skip empty lines

        try:
            # 1. Generate RDKit Mol
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                print(f'    Warning: Invalid SMILES at line {i+1} in {base_name}')
                continue
            
            # 2. Add Hydrogens & Embed 3D
            mol = Chem.AddHs(mol)
            res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            if res != 0:
                print(f'    Warning: 3D embedding failed for line {i+1} in {base_name}')
                continue

            # 3. Prepare PDBQT with Meeko
            prep = MoleculePreparation()
            setups = prep.prepare(mol)
            
            data = PDBQTWriterLegacy.write_string(setups[0])
            string = data[0] if isinstance(data, tuple) else data
            
            # 4. Write to unique file (Glutamine_0.pdbqt, Glutamine_1.pdbqt...)
            out_filename = f'{base_name}_{count}.pdbqt'
            out_path = os.path.join(out_dir, out_filename)
            
            with open(out_path, 'w') as f:
                f.write(string)
            
            count += 1

        except Exception as inner_e:
            print(f'    Error processing molecule {i} in {base_name}: {inner_e}')

    print(f'    > Generated {count} PDBQT files for {base_name}')

except Exception as e:
    print(f'Error reading file {sml_file}: {e}')
" "$LIGAND_SML" "$LIGAND_PDBQT_DIR" "$LIG_NAME"

done

# Initialize Summary File
SUMMARY_FILE="$OUTPUT_BASE_DIR/summary_scores.csv"
echo "Receptor,Pocket_ID,Ligand,Affinity" > "$SUMMARY_FILE"


# --- PHASE 2: TARGETED DOCKING (BioBB / fpocket) ---
echo "=== Phase 2: Targeted Docking (BioBB/fpocket) ==="

for RECEPTOR_PDB in "$STRUCTURES_DIR"/*.pdb; do
    REC_NAME=$(basename -- "$RECEPTOR_PDB" .pdb)
    RECEPTOR_PDBQT="$STRUCTURES_DIR/${REC_NAME}.pdbqt"
    
    # Temp folder for BioBB calculations (cleaned up later if desired)
    BIOBB_TEMP="$OUTPUT_BASE_DIR/${REC_NAME}_biobb"

    # 1. Prep Receptor PDBQT
    if [ ! -f "$RECEPTOR_PDBQT" ]; then
        obabel -ipdb "$RECEPTOR_PDB" -opdbqt -O "$RECEPTOR_PDBQT" -xr -h --partialcharge gasteiger
    fi

    # 2. Find Pockets
    # Calls the external python script to get pocket centers
    POCKET_DATA=$(python3 "$BASE_DIR/scripts/BioBB_pocket_finder.py" "$RECEPTOR_PDB" "$BIOBB_TEMP")

    if [ -z "$POCKET_DATA" ]; then
        echo "  > No pockets found for $REC_NAME matching criteria."
        continue
    fi

    IFS=$'\n'
    for POCKET_LINE in $POCKET_DATA; do
        # 1. Clean the line (remove extra whitespace)
        CLEAN_LINE=$(echo "$POCKET_LINE" | xargs)
        
        # 2. Skip if empty
        if [ -z "$CLEAN_LINE" ]; then continue; fi

        # 3. Read variables
        read -r P_ID P_SCORE CX CY CZ <<< "$CLEAN_LINE"
        
        # 4. SAFETY CHECK: Ensure P_ID is actually a number.
        # If BioBB printed an error message, P_ID might be "Error" or "Warning"
        if ! [[ "$P_ID" =~ ^[0-9]+$ ]]; then
            echo "  > Skipping non-data line from pocket finder: $CLEAN_LINE"
            continue
        fi
        
        # Create output folder: outputs/ReceptorName/pocket_1/
        OUT_DIR="$OUTPUT_BASE_DIR/${REC_NAME}/pocket_${P_ID}"
        mkdir -p "$OUT_DIR"
        
        echo "  > Docking into ${REC_NAME} Pocket $P_ID (Center: $CX $CY $CZ)"

        # 3. Dock All Ligands into this Pocket
        for LIGAND_PDBQT in "$LIGAND_PDBQT_DIR"/*.pdbqt; do
            LIG_NAME=$(basename -- "$LIGAND_PDBQT" .pdbqt)
            OUT_POSE="$OUT_DIR/${LIG_NAME}_docked.pdbqt"
            LOG="$OUT_DIR/${LIG_NAME}.log"

            vina --receptor "$RECEPTOR_PDBQT" --ligand "$LIGAND_PDBQT" \
                 --center_x $CX --center_y $CY --center_z $CZ \
                 --size_x $SIZE_X --size_y $SIZE_Y --size_z $SIZE_Z \
                 --exhaustiveness $EXHAUSTIVENESS \
                 --out "$OUT_POSE" > "$LOG" 2>&1

            # 4. Parse Score
            if [ -f "$LOG" ]; then
                SCORE=$(awk '/^ *1 / {print $2; exit}' "$LOG")
                # Append to summary CSV
                echo "$REC_NAME,$P_ID,$LIG_NAME,$SCORE" >> "$SUMMARY_FILE"
            fi
        done
    done
    unset IFS
done

echo "========================================"
echo "Targeted Docking Complete."
echo "Summary saved to: $SUMMARY_FILE"