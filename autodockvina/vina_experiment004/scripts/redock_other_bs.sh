#!/bin/bash

# ==============================================================================
# VINA SPECIFICITY CHECK PROTOCOL
# ==============================================================================

# --- CONFIGURATION ---
BASE_DIR="/workspace/vina_experiment004"

# 1. RECEPTOR SELECTION
# List specific filenames (without .pdb extension).
# These are used to filter the loops and locate the specific Pockets CSV.
SELECTED_RECEPTORS=("SLC6A20_8I91_complete_occluded")

# 2. INPUT DATA
# We use the first element of the array to define the input filenames
# (Assuming the input CSV is named after the main receptor target)
MAIN_TARGET="${SELECTED_RECEPTORS[0]}"

# The Pockets CSV (Output from p2rank/fpocket comparison)
# MUST have columns: 'Receptor', 'p2rank_rank', 'fpocket_nr', 'p2Rank_center', 'fpocket_center'
POCKETS_CSV="$BASE_DIR/data/pockets/${MAIN_TARGET}_pockets_final.csv"

# The results from the specificity ranking script
# Contains: Unique_ID, Rank_A, Rank_B, Affinity_A, Affinity_B, Rank_Delta
NATURAL_RESULTS_CSV="$BASE_DIR/outputs/orthosteric_auto/receptor_comparison/specificity_SLC6A20__occluded__vs_SLC6A19__occluded__hits.csv"

# Paths
STRUCTURES_DIR="$BASE_DIR/data/structures"
LIGAND_LIBRARY_DIR="$BASE_DIR/data/library1_pdbqt"

# 3. OUTPUT
OUTPUT_DIR="$BASE_DIR/outputs/redocking"
SUMMARY_FILE="$OUTPUT_DIR/specificity_scores.csv"

# 4. SETTINGS
TOP_N=10         # Number of top ligands to test
GRID_SIZE=20     # Box size for new pockets

# --- SETUP ---
mkdir -p "$OUTPUT_DIR"

# Write Header if file doesn't exist
if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Receptor,Pocket_ID,Ligand_ID,Affinity_kcal_mol" > "$SUMMARY_FILE"
fi

echo "========================================"
echo "SPECIFICITY CHECK PROTOCOL"
echo "Target Group: $MAIN_TARGET"
echo "Receptors to process: ${SELECTED_RECEPTORS[*]}"
echo "========================================"

# --- STEP 1: GET TOP HITS (Once) ---
echo "1. Extracting Top $TOP_N Specific Ligands..."

# Export variables for Python to read safely
export NATURAL_RESULTS_CSV
export TOP_N

python3 -c "
import pandas as pd
import sys
import os

csv_path = os.environ['NATURAL_RESULTS_CSV']
top_n = int(os.environ['TOP_N'])

try:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f'File not found: {csv_path}')

    df = pd.read_csv(csv_path)
    
    # Sort by Rank_Delta if available (highest/best on top)
    if 'Rank_Delta' in df.columns:
        df_sorted = df.sort_values(by='Rank_Delta', ascending=False)
    else:
        df_sorted = df

    top_hits = df_sorted.head(top_n)
    
    # Identify ID column
    if 'Unique_ID' in top_hits.columns:
        id_col = 'Unique_ID'
    elif 'Ligand_ID' in top_hits.columns:
        id_col = 'Ligand_ID'
    else:
         raise ValueError('Unique_ID column not found in CSV')

    for val in top_hits[id_col]:
        print(val)

except Exception as e:
    print(f'Error reading CSV: {e}', file=sys.stderr)
    sys.exit(1)
" > "$OUTPUT_DIR/top_ligands.txt"

# Check if Python succeeded
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract top ligands."
    exit 1
fi

mapfile -t TOP_LIGANDS < "$OUTPUT_DIR/top_ligands.txt"

if [ ${#TOP_LIGANDS[@]} -eq 0 ]; then
    echo "Error: No top ligands found. Check CSV path: $NATURAL_RESULTS_CSV"
    exit 1
fi
echo "   -> Selected ligands: ${TOP_LIGANDS[*]}"


# --- STEP 2: LOOP OVER RECEPTORS ---
# We loop through all PDBs in the directory, but only process those in the filtered list
for RECEPTOR_PDB in "$STRUCTURES_DIR"/*.pdb; do
    
    REC_FILENAME=$(basename -- "$RECEPTOR_PDB")
    REC_NAME="${REC_FILENAME%.*}"
    RECEPTOR_PDBQT="$STRUCTURES_DIR/${REC_NAME}.pdbqt"

    # --- SELECTION LOGIC ---
    PROCESS_THIS=false
    if [ "${SELECTED_RECEPTORS[0]}" == "all" ]; then
        PROCESS_THIS=true
    else
        for TARGET in "${SELECTED_RECEPTORS[@]}"; do
            if [[ "$REC_NAME" == "$TARGET" ]]; then
                PROCESS_THIS=true; break
            fi
        done
    fi

    if [ "$PROCESS_THIS" = false ]; then continue; fi
    # -----------------------

    echo "========================================"
    echo "PROCESSING RECEPTOR: $REC_NAME"
    
    # Convert PDB to PDBQT if needed
    if [ ! -f "$RECEPTOR_PDBQT" ]; then
        echo "   -> Converting PDB to PDBQT..."
        obabel -ipdb "$RECEPTOR_PDB" -opdbqt -O "$RECEPTOR_PDBQT" -xr -h --partialcharge gasteiger 2> /dev/null
        if [ ! -f "$RECEPTOR_PDBQT" ]; then
            echo "Error: Obabel conversion failed for $REC_NAME"
            continue
        fi
    fi

    # Create Receptor Output Folder
    REC_OUT_DIR="$OUTPUT_DIR/$REC_NAME"
    mkdir -p "$REC_OUT_DIR"

    # --- STEP 3: LOOP OVER POCKETS (Filtered by Receptor) ---
    export POCKETS_CSV
    export REC_NAME
    
    # We use a temporary file for pocket data to avoid subshell variable loss issues
    POCKET_DATA_FILE="$REC_OUT_DIR/pockets_to_dock.txt"

    python3 -c "
import pandas as pd
import sys
import re
import os

csv_path = os.environ['POCKETS_CSV']
target_rec = os.environ['REC_NAME']

def get_coords(row, col_name):
    val = str(row.get(col_name, ''))
    # Remove brackets [] () and split by comma or whitespace
    # Using raw string for regex to be safe
    cleaned = re.sub(r'[\[\]\(\)]', '', val)
    parts = re.split(r'[,\s]+', cleaned.strip())
    parts = [p for p in parts if p]
    
    if len(parts) >= 3:
        try:
            return float(parts[0]), float(parts[1]), float(parts[2])
        except:
            return None
    return None

try:
    if not os.path.exists(csv_path):
        # Fail silently here if file missing, handled by bash
        sys.exit(0)

    df = pd.read_csv(csv_path)
    
    # Iterate rows and filter by Receptor
    for index, row in df.iterrows():
        csv_rec = str(row.get('Receptor', '')).strip()
        
        # Check matching (partial match logic)
        if csv_rec == target_rec or csv_rec in target_rec or target_rec in csv_rec:
            
            p2_rank = row.get('p2rank_rank', 'X')
            fp_nr = row.get('fpocket_nr', 'X')
            pid = f'p2r{p2_rank}_fp{fp_nr}'
            
            # Extract Coordinates
            xyz = get_coords(row, 'p2Rank_center')
            if not xyz:
                xyz = get_coords(row, 'fpocket_center')
            
            if xyz:
                # Print: ID X Y Z
                print(f'{pid} {xyz[0]} {xyz[1]} {xyz[2]}')

except Exception as e:
    print(f'Error parsing pockets: {e}', file=sys.stderr)
" > "$POCKET_DATA_FILE"

    if [ ! -s "$POCKET_DATA_FILE" ]; then
        echo "   > No pockets found for $REC_NAME (or CSV missing)."
        continue
    fi

    # Process the pockets found
    while read -r POCKET_ID CX CY CZ; do

        # Sanity check coordinates
        if [[ -z "$CX" || -z "$CY" || -z "$CZ" ]]; then
            echo "   > Skipping invalid coords for $POCKET_ID"
            continue
        fi

        echo "   > Pocket: $POCKET_ID (Center: $CX $CY $CZ)"
        
        POCKET_DIR="$REC_OUT_DIR/$POCKET_ID"
        mkdir -p "$POCKET_DIR"

        # --- STEP 4: DOCK TOP LIGANDS ---
        for LIG_RAW in "${TOP_LIGANDS[@]}"; do
            
            LIG_ID="$LIG_RAW"
            LIG_PATH=""

            # HEURISTIC FILE FINDING
            TEMP_ID="$LIG_ID"
            
            # 1. Exact match
            if [ -f "$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt" ]; then
                LIG_PATH="$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt"
            else
                # 2. Strip prefix iteratively (e.g. Lib_Index -> Index)
                while [[ "$TEMP_ID" == *"_"* ]]; do
                    TEMP_ID="${TEMP_ID#*_}"
                    if [ -f "$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt" ]; then
                        LIG_PATH="$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt"
                        break
                    fi
                done
            fi
            
            # 3. Find by suffix
            if [ -z "$LIG_PATH" ]; then
                 CLEAN_ID=$(echo "$LIG_ID" | awk -F_ '{print $NF}')
                 LIG_PATH=$(find "$LIGAND_LIBRARY_DIR" -name "*${CLEAN_ID}.pdbqt" | head -n 1)
            fi

            if [ -z "$LIG_PATH" ]; then
                echo "      Warning: PDBQT file for $LIG_ID not found."
                continue
            fi

            LIG_BASENAME=$(basename "$LIG_PATH" .pdbqt)
            OUT_POSE="$POCKET_DIR/${LIG_BASENAME}_redock.pdbqt"
            LOG_FILE="$POCKET_DIR/${LIG_BASENAME}_redock.log"
            
            # Check Summary to skip done work
            if grep -q "$REC_NAME,$POCKET_ID,$LIG_BASENAME" "$SUMMARY_FILE"; then
                echo "      - Skipping $LIG_BASENAME (already in summary)"
                continue
            fi

            # Docking
            vina \
                --receptor "$RECEPTOR_PDBQT" \
                --ligand "$LIG_PATH" \
                --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
                --size_x "$GRID_SIZE" --size_y "$GRID_SIZE" --size_z "$GRID_SIZE" \
                --exhaustiveness 8 \
                --out "$OUT_POSE" \
                > "$LOG_FILE" 2>&1

            if [ -f "$LOG_FILE" ]; then
                # Robust Grep for Vina Score (Mode 1)
                # Looks for line starting with "   1" followed by numbers
                SCORE=$(grep "^ \+1 " "$LOG_FILE" | awk '{print $2}')
                
                if [ -n "$SCORE" ]; then
                    echo "$REC_NAME,$POCKET_ID,$LIG_BASENAME,$SCORE" >> "$SUMMARY_FILE"
                else
                    echo "      ! Docking failed or no score found for $LIG_BASENAME"
                fi
            fi

        done # End Ligands
    done < "$POCKET_DATA_FILE" # End Pockets loop

done # End Receptors

echo "Specificity Check Complete. Results in $SUMMARY_FILE"