#!/bin/bash

# --- CONFIGURATION ---
BASE_DIR="/workspace/vina_experiment002"

# What and where should it dock?
TARGET_SITE="orthosteric"   # Optionen: "orthosteric", "allosteric", "all"
TARGET_CONF="auto"          # Optionen: "occluded", "inward", "outward", "auto" (nimmt was im Dateinamen steht)

# Binding pocket Residues 
BINDING_DATA_FILE="$BASE_DIR/data/reference/pockets.csv"

# Input Directories
STRUCTURES_DIR="$BASE_DIR/data/structures"
LIGAND_SML_DIR="$BASE_DIR/data/library1_smiles"

# Intermediate/Output Directories
LIGAND_PDBQT_DIR="$BASE_DIR/data/library1_pdbqt"
OUTPUT_BASE_DIR="$BASE_DIR/outputs/${TARGET_SITE}_${TARGET_CONF}"

# Vina Grid Size (Angstroms)
GRID_SIZE=20
EXHAUSTIVENESS=8

# --- SETUP ---
mkdir -p "$OUTPUT_BASE_DIR"
PYTHON_SCRIPT="$BASE_DIR/scripts/get_pocket_center.py"
SUMMARY_FILE="$OUTPUT_BASE_DIR/summary_scores.csv"

if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Receptor,Ligand_Library,Ligand_ID,Affinity_kcal_mol" > "$SUMMARY_FILE"
fi

echo "========================================"
echo "STARTING DOCKING PIPELINE"
echo "Target Site: $TARGET_SITE"
echo "Target Conf: $TARGET_CONF"
echo "========================================"

# --- OUTER LOOP: RECEPTORS ---
for RECEPTOR_PDB in "$STRUCTURES_DIR"/*.pdb; do
    
    REC_FILENAME=$(basename -- "$RECEPTOR_PDB")
    REC_NAME="${REC_FILENAME%.*}"
    RECEPTOR_PDBQT="$STRUCTURES_DIR/${REC_NAME}.pdbqt"

    # 1. CHECK AVAILABILITY & CALCULATE CENTER
    # Wir rufen Python auf mit den User-Wünschen
    # Wenn Python nichts zurückgibt (leer), überspringen wir den Rezeptor
    
    COORDS=$(python3 "$PYTHON_SCRIPT" "$RECEPTOR_PDB" "$BINDING_DATA_FILE" "$REC_NAME" --site "$TARGET_SITE" --conf "$TARGET_CONF")

    if [ -z "$COORDS" ]; then
        echo "  > Skipping $REC_NAME: No binding site found for $TARGET_SITE in $TARGET_CONF conformation."
        continue
    fi

    # Koordinaten lesen
    read CENTER_X CENTER_Y CENTER_Z <<< "$COORDS"
    
    echo "  > Processing $REC_NAME"
    echo "    Pocket found at: $CENTER_X $CENTER_Y $CENTER_Z"

    # Ordner erstellen
    REC_OUTPUT_DIR="$OUTPUT_BASE_DIR/$REC_NAME"
    mkdir -p "$REC_OUTPUT_DIR"

    # PDBQT Konvertierung (falls nötig)
    if [ ! -f "$RECEPTOR_PDBQT" ]; then
        obabel -ipdb "$RECEPTOR_PDB" -opdbqt -O "$RECEPTOR_PDBQT" -xr -h --partialcharge gasteiger 2> /dev/null
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
        out_filename = f'{base_name}_{idx:04d}.pdbqt'
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
                --size_x $GRID_SIZE --size_y $GRID_SIZE --size_z $GRID_SIZE \
                --exhaustiveness $EXHAUSTIVENESS \
                --out "$OUTPUT_POSE" \
                > "$LOG_FILE" 2>&1

            # Parse Score
            if [ -f "$LOG_FILE" ]; then
                AFFINITY=$(grep "   1 " "$LOG_FILE" | head -n 1 | awk '{print $2}')
                # Schreibt jetzt in die korrekt definierte Datei
                echo "$REC_NAME,$LIG_NAME,$SUB_LIG_ID,$AFFINITY" >> "$SUMMARY_FILE"
            fi
        
        done # End Sub-Ligand Loop

    done # End SML File Loop

done # End Receptor Loop
echo "Batch Docking Complete."