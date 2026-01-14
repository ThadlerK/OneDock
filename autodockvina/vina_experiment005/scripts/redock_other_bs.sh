#!/bin/bash

# --- CONFIGURATION ---
BASE_DIR="/workspace/vina_experiment005"

# 1. RECEPTOR SELECTION
# Set to ("all") to dock against every PDB in the structures directory.
# OR list specific filenames (without .pdb extension) to limit the run.
# IMPORTANT: These names must match (or be contained in) the "Receptor" column of your pockets CSV.
SELECTED_RECEPTORS=("SLC6A20_8WM3_complete_inwardopen")

# 2. INPUT DATA
# The Pockets CSV (Output from p2rank/fpocket comparison)
# MUST have columns: 'Receptor', 'p2rank_rank', 'fpocket_nr', 'p2Rank_center', 'fpocket_center'
POCKETS_CSV="$BASE_DIR/data/pockets/pockets_final.csv"

# The results from the specificity ranking script
# Contains: Unique_ID, Rank_A, Rank_B, Affinity_A, Affinity_B, Rank_Delta
NATURAL_RESULTS_CSV="$BASE_DIR/outputs/analysis_ranking/specific_hits.csv"

# Paths
STRUCTURES_DIR="$BASE_DIR/data/structures"
LIGAND_LIBRARY_DIR="$BASE_DIR/data/library1_pdbqt"

# 3. OUTPUT
OUTPUT_DIR="$BASE_DIR/outputs/specificity_check"
SUMMARY_FILE="$OUTPUT_DIR/specificity_scores.csv"

# 4. SETTINGS
TOP_N=10         # Number of top ligands to test
GRID_SIZE=20     # Box size for new pockets

# --- SETUP ---
mkdir -p "$OUTPUT_DIR"

# Write Header (including Receptor column)
if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Receptor,Pocket_ID,Ligand_ID,Affinity_kcal_mol" > "$SUMMARY_FILE"
fi

echo "========================================"
echo "SPECIFICITY CHECK PROTOCOL"
echo "Receptors: ${SELECTED_RECEPTORS[*]}"
echo "========================================"

# --- STEP 1: GET TOP HITS (Once) ---
echo "1. Extracting Top $TOP_N Specific Ligands..."

python3 -c "
import pandas as pd
import sys

try:
    df = pd.read_csv('$NATURAL_RESULTS_CSV')
    
    # User states file is already sorted by Rank_Delta (highest/best on top).
    # We take the top N rows directly.
    # Optionally, we can enforce sort if Rank_Delta exists.
    if 'Rank_Delta' in df.columns:
        df_sorted = df.sort_values(by='Rank_Delta', ascending=False)
    else:
        # Fallback: assume file is already sorted
        df_sorted = df

    top_hits = df_sorted.head($TOP_N)
    
    # Identify ID column (Unique_ID is standard from ranking script)
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

mapfile -t TOP_LIGANDS < "$OUTPUT_DIR/top_ligands.txt"

if [ ${#TOP_LIGANDS[@]} -eq 0 ]; then
    echo "Error: No top ligands found. Check CSV path!"
    exit 1
fi
echo "   -> Selected ligands: ${TOP_LIGANDS[*]}"


# --- STEP 2: LOOP OVER RECEPTORS ---
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
    
    # Convert if needed
    if [ ! -f "$RECEPTOR_PDBQT" ]; then
        obabel -ipdb "$RECEPTOR_PDB" -opdbqt -O "$RECEPTOR_PDBQT" -xr -h --partialcharge gasteiger 2> /dev/null
    fi

    # Create Receptor Output Folder
    REC_OUT_DIR="$OUTPUT_DIR/$REC_NAME"
    mkdir -p "$REC_OUT_DIR"

    # --- STEP 3: LOOP OVER POCKETS (Filtered by Receptor) ---
    # We use Python to parse the CSV and only return pockets 
    # where the 'Receptor' column matches the current REC_NAME.
    
    POCKET_DATA=$(python3 -c "
import pandas as pd
import sys
import re

def get_coords(row, col_name):
    val = str(row.get(col_name, ''))
    # Remove brackets [] () and split by comma or whitespace
    cleaned = re.sub(r'[\[\]\(\)]', '', val)
    parts = re.split(r'[,\s]+', cleaned.strip())
    # Filter empty strings
    parts = [p for p in parts if p]
    
    if len(parts) >= 3:
        try:
            return float(parts[0]), float(parts[1]), float(parts[2])
        except:
            return None
    return None

try:
    df = pd.read_csv('$POCKETS_CSV')
    
    target_rec = '$REC_NAME'
    found_any = False

    # Iterate rows and filter by Receptor
    for index, row in df.iterrows():
        # Clean up CSV entry (remove whitespace)
        csv_rec = str(row.get('Receptor', '')).strip()
        
        # Check if names match (allows partial matches like 'SLC6A20' in 'SLC6A20_occluded')
        if csv_rec == target_rec or csv_rec in target_rec or target_rec in csv_rec:
            
            # Construct a clean Pocket ID
            p2_rank = row.get('p2rank_rank', 'X')
            fp_nr = row.get('fpocket_nr', 'X')
            pid = f'p2r{p2_rank}_fp{fp_nr}'
            
            # Extract Coordinates (Prefer p2Rank, fallback to fpocket)
            xyz = get_coords(row, 'p2Rank_center')
            if not xyz:
                xyz = get_coords(row, 'fpocket_center')
            
            if xyz:
                # Print: ID X Y Z
                print(f'{pid} {xyz[0]} {xyz[1]} {xyz[2]}')
                found_any = True

except Exception as e:
    print(f'Error parsing pockets: {e}', file=sys.stderr)
")

    # If no pockets found for this receptor, skip docking
    if [ -z "$POCKET_DATA" ]; then
        echo "   > No pockets found for $REC_NAME in $POCKETS_CSV. Skipping."
        continue
    fi

    # Process the pockets found
    echo "$POCKET_DATA" | while read POCKET_ID CX CY CZ; do

        echo "   > Pocket: $POCKET_ID"
        
        POCKET_DIR="$REC_OUT_DIR/$POCKET_ID"
        mkdir -p "$POCKET_DIR"

        # --- STEP 4: DOCK TOP LIGANDS ---
        for LIG_RAW in "${TOP_LIGANDS[@]}"; do
            
            LIG_ID="$LIG_RAW"
            LIG_PATH=""

            # HEURISTIC FILE FINDING
            # Unique_ID is often "Library_Library_Index" while filename is "Library_Index.pdbqt"
            # We strip prefixes iteratively until we find a matching file.
            
            TEMP_ID="$LIG_ID"
            # 1. Check direct match first
            if [ -f "$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt" ]; then
                LIG_PATH="$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt"
            else
                # 2. Iteratively strip prefix (everything up to first underscore)
                while [[ "$TEMP_ID" == *"_"* ]]; do
                    TEMP_ID="${TEMP_ID#*_}"
                    if [ -f "$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt" ]; then
                        LIG_PATH="$LIGAND_LIBRARY_DIR/${TEMP_ID}.pdbqt"
                        break
                    fi
                done
            fi
            
            # 3. Last resort fallback (find by suffix)
            if [ -z "$LIG_PATH" ]; then
                 CLEAN_ID=$(echo "$LIG_ID" | awk -F_ '{print $NF}')
                 LIG_PATH=$(find "$LIGAND_LIBRARY_DIR" -name "*${CLEAN_ID}.pdbqt" | head -n 1)
            fi

            if [ -z "$LIG_PATH" ]; then
                echo "      Warning: PDBQT file for $LIG_ID not found."
                continue
            fi

            # Clean Name for Output Files (e.g. Library_0001)
            LIG_BASENAME=$(basename "$LIG_PATH" .pdbqt)

            OUT_POSE="$POCKET_DIR/${LIG_BASENAME}_redock.pdbqt"
            LOG_FILE="$POCKET_DIR/${LIG_BASENAME}_redock.log"
            
            # Check Summary (avoids double work)
            if grep -q "$REC_NAME,$POCKET_ID,$LIG_BASENAME" "$SUMMARY_FILE"; then
                continue
            fi

            # Docking
            vina \
                --receptor "$RECEPTOR_PDBQT" \
                --ligand "$LIG_PATH" \
                --center_x $CX --center_y $CY --center_z $CZ \
                --size_x $GRID_SIZE --size_y $GRID_SIZE --size_z $GRID_SIZE \
                --exhaustiveness 8 \
                --out "$OUT_POSE" \
                > "$LOG_FILE" 2>&1

            if [ -f "$LOG_FILE" ]; then
                SCORE=$(grep "   1 " "$LOG_FILE" | head -n 1 | awk '{print $2}')
                # Save: Receptor, Pocket, Ligand, Score
                echo "$REC_NAME,$POCKET_ID,$LIG_BASENAME,$SCORE" >> "$SUMMARY_FILE"
            fi

        done # End Ligands
    done # End Pockets
done # End Receptors

echo "Specificity Check Complete. Results in $SUMMARY_FILE"