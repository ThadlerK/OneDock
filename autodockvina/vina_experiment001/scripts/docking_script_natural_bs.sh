#!/bin/bash

# --- CONFIGURATION ---
BASE_DIR="/workspace/vina_experiment002"

# Reference strucutre with ligand bound
REFERENCE_PDB="/workspace/vina_experiment002/data/reference/8I91.pdb"
# Ligand name
REF_LIGAND_NAME="PRO"

# Input Directories
STRUCTURES_DIR="$BASE_DIR/data/structures"
LIGAND_SML_DIR="$BASE_DIR/data/library1_smiles"

# Intermediate/Output Directories
LIGAND_PDBQT_DIR="$BASE_DIR/data/library1_pdbqt"
OUTPUT_BASE_DIR="$BASE_DIR/outputs"

# Vina Grid Size (Angstroms)
GRID_SIZE_X=18
GRID_SIZE_Y=18
GRID_SIZE_Z=18
EXHAUSTIVENESS=8

# --- SETUP ---
mkdir -p "$LIGAND_PDBQT_DIR"
mkdir -p "$OUTPUT_BASE_DIR"
SUMMARY_FILE="$OUTPUT_BASE_DIR/summary_scores.csv"
if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Receptor,Ligand_File,Ligand_Index,Affinity_kcal_mol" > "$SUMMARY_FILE"
fi

# --- STEP 1: CALCULATE GRID CENTER (GLOBAL) ---
echo "========================================"
echo "DEFINING POCKET BASED ON: $REFERENCE_PDB ($REF_LIGAND_NAME)"

read CENTER_X CENTER_Y CENTER_Z <<< $(grep "$REF_LIGAND_NAME" "$REFERENCE_PDB" | awk '$1=="ATOM" || $1=="HETATM" {x+=$7; y+=$8; z+=$9; n++} END {printf "%.3f %.3f %.3f", x/n, y/n, z/n}')

if [ -z "$CENTER_X" ]; then
    echo "Error: Could not find ligand $REF_LIGAND_NAME in $REFERENCE_PDB"
    exit 1
fi
echo "  > Pocket Center: $CENTER_X, $CENTER_Y, $CENTER_Z"
echo "========================================"

# --- OUTER LOOP: RECEPTORS ---
for RECEPTOR_PDB in "$STRUCTURES_DIR"/*.pdb; do
    
    REC_FILENAME=$(basename -- "$RECEPTOR_PDB")
    REC_NAME="${REC_FILENAME%.*}"
    RECEPTOR_PDBQT="$STRUCTURES_DIR/${REC_NAME}.pdbqt"
    REC_OUTPUT_DIR="$OUTPUT_BASE_DIR/$REC_NAME"
    mkdir -p "$REC_OUTPUT_DIR"

    echo "========================================"
    echo "PROCESSING RECEPTOR: $REC_NAME"

    # Prepare Receptor
    if [ ! -f "$RECEPTOR_PDBQT" ]; then
        echo "  > Converting Receptor to PDBQT..."
        obabel -ipdb "$RECEPTOR_PDB" -opdbqt -O "$RECEPTOR_PDBQT" -xr -h --partialcharge gasteiger
    fi

    # --- MIDDLE LOOP: SML FILES ---
    for LIGAND_SML in "$LIGAND_SML_DIR"/*.sml; do
        
        LIG_FILENAME=$(basename -- "$LIGAND_SML")
        LIG_NAME="${LIG_FILENAME%.*}" # z.B. "Glutamine"

        echo "    ------------------------------------"
        echo "    Processing SML File: $LIG_NAME"

        # 1. PYTHON: SPLIT & CONVERT
        # Wir übergeben jetzt das Output-Verzeichnis, nicht einen einzelnen Dateinamen
        # Python generiert: LIG_NAME_0.pdbqt, LIG_NAME_1.pdbqt, etc.
        
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
        # Lese alle Zeilen, ignoriere leere Zeilen
        lines = [l.strip() for l in f if l.strip()]

    print(f'      > Found {len(lines)} SMILES in file.')

    for idx, line in enumerate(lines):
        # Nimm nur den ersten Teil (falls Namen hinter dem SMILES stehen)
        smiles = line.split()[0]
        
        # Output Name definieren: Glutamine_0.pdbqt, Glutamine_1.pdbqt...
        out_filename = f'{base_name}_{idx}.pdbqt'
        out_path = os.path.join(out_dir, out_filename)

        # Check ob schon existiert (um Zeit zu sparen)
        if os.path.exists(out_path):
            continue

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f'      > Warning: Invalid SMILES at line {idx+1}')
            continue

        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            pdbqt_data = PDBQTWriterLegacy.write_string(mol_setups[0])
            pdbqt_string = pdbqt_data[0] if isinstance(pdbqt_data, tuple) else pdbqt_data

            with open(out_path, 'w') as out_f:
                out_f.write(pdbqt_string)
        except Exception as e:
            print(f'      > Error prepping ligand {idx}: {e}')

except Exception as e:
    print(f'Error opening {sml_file}: {e}')
    sys.exit(1)
" "$LIGAND_SML" "$LIGAND_PDBQT_DIR" "$LIG_NAME"


        # --- INNER LOOP: DOCK GENERATED PDBQT FILES ---
        # Jetzt suchen wir alle Dateien, die Python gerade erstellt hat (Glutamine_*.pdbqt)
        
        # nullglob verhindert Fehler, falls keine Dateien gefunden werden
        shopt -s nullglob
        PDBQT_FILES=("$LIGAND_PDBQT_DIR/${LIG_NAME}_"*.pdbqt)
        shopt -u nullglob

        if [ ${#PDBQT_FILES[@]} -eq 0 ]; then
            echo "      > No valid PDBQT files generated for $LIG_NAME"
            continue
        fi

        for SUB_LIGAND_PATH in "${PDBQT_FILES[@]}"; do
            
            # Filename extrahieren: "Glutamine_0.pdbqt"
            SUB_FILENAME=$(basename -- "$SUB_LIGAND_PATH")
            # Name ohne Endung: "Glutamine_0"
            SUB_LIG_ID="${SUB_FILENAME%.*}" 
            
            OUTPUT_POSE="$REC_OUTPUT_DIR/${SUB_LIG_ID}_docked.pdbqt"
            LOG_FILE="$REC_OUTPUT_DIR/${SUB_LIG_ID}_docked.log"

            # Check if already docked
            if [ -f "$SUMMARY_FILE" ] && grep -q "$REC_NAME,$LIG_NAME,$SUB_LIG_ID" "$SUMMARY_FILE"; then
                 # echo "Skipping $SUB_LIG_ID (already done)"
                 continue
            fi

            echo "      > Docking variant: $SUB_LIG_ID"

            vina \
                --receptor "$RECEPTOR_PDBQT" \
                --ligand "$SUB_LIGAND_PATH" \
                --center_x $CENTER_X \
                --center_y $CENTER_Y \
                --center_z $CENTER_Z \
                --size_x $GRID_SIZE_X --size_y $GRID_SIZE_Y --size_z $GRID_SIZE_Z \
                --exhaustiveness $EXHAUSTIVENESS \
                --out "$OUTPUT_POSE" \
                > "$LOG_FILE" 2>&1

            # Parse Score
            if [ -f "$LOG_FILE" ]; then
                AFFINITY=$(grep "   1 " "$LOG_FILE" | head -n 1 | awk '{print $2}')
                # CSV Format: Receptor, SML_File_Name, Specific_Ligand_ID, Score
                echo "$REC_NAME,$LIG_NAME,$SUB_LIG_ID,$AFFINITY" >> "$SUMMARY_FILE"
            fi
        
        done # End Sub-Ligand Loop

    done # End SML File Loop

done # End Receptor Loop

echo "Batch Docking Complete."