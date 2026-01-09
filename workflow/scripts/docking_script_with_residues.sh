#!/bin/bash

# --- INPUT ARGUMENTS ---
RECEPTOR_PDBQT="$1"
LIGAND_PDBQT="$2"
RESIDUE_LIST="$3"       # Expected format: "A:123,A:124" or just "123,124" if no chain
OUTPUT_DOCKED="$4"
OUTPUT_LOG="$5"
OUTPUT_SUMMARY="$6"
GRID_SIZE="${7:-20}"    # Default to 20 if not set
EXHAUSTIVENESS="${8:-8}" # Default to 8 if not set

echo "========================================"
echo "STARTING SINGLE DOCKING RUN"
echo "Receptor: $RECEPTOR_PDBQT"
echo "Ligand:   $LIGAND_PDBQT"
echo "Residues: $RESIDUE_LIST"
echo "========================================"

# Create directories
mkdir -p "$(dirname "$OUTPUT_DOCKED")"
mkdir -p "$(dirname "$OUTPUT_LOG")"
mkdir -p "$(dirname "$OUTPUT_SUMMARY")"

# --- 1. CALCULATE POCKET CENTER FROM RESIDUES ---
# We use a quick inline Python script to parse the PDBQT and find the center of mass of the selected residues.
# This ensures the center is calculated exactly for the input file used.

COORDS=$(python3 -c "
import sys

pdbqt_file = sys.argv[1]
residue_str = sys.argv[2] # Format 'A:123,B:456' or '123,456'

# Parse target residues
targets = set()
for item in residue_str.split(','):
    if ':' in item:
        chain, resnum = item.split(':')
        targets.add((chain, int(resnum)))
    else:
        # Handle case where no chain is provided (assume matching any chain or blank)
        targets.add(('*', int(item)))

coords = []

with open(pdbqt_file, 'r') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                # PDBQT column positions (1-based index):
                # ChainID: col 22 (index 21)
                # ResSeq:  col 23-26 (index 22-26)
                # X: 31-38, Y: 39-46, Z: 47-54
                
                chain_id = line[21].strip()
                res_seq = int(line[22:26].strip())
                
                match = False
                if (chain_id, res_seq) in targets:
                    match = True
                elif ('*', res_seq) in targets:
                    match = True
                
                if match:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
            except ValueError:
                continue

if not coords:
    print('ERROR')
else:
    # Calculate geometric center
    avg_x = sum(c[0] for c in coords) / len(coords)
    avg_y = sum(c[1] for c in coords) / len(coords)
    avg_z = sum(c[2] for c in coords) / len(coords)
    print(f'{avg_x:.3f} {avg_y:.3f} {avg_z:.3f}')
" "$RECEPTOR_PDBQT" "$RESIDUE_LIST")

if [ "$COORDS" == "ERROR" ] || [ -z "$COORDS" ]; then
    echo "Error: Could not calculate center. Check if residues $RESIDUE_LIST exist in $RECEPTOR_PDBQT"
    exit 1
fi

read CENTER_X CENTER_Y CENTER_Z <<< "$COORDS"
echo "Calculated Center: $CENTER_X $CENTER_Y $CENTER_Z"

# --- 2. RUN VINA ---
# Ensure output directory exists
mkdir -p "$(dirname "$OUTPUT_DOCKED")"

vina \
    --receptor "$RECEPTOR_PDBQT" \
    --ligand "$LIGAND_PDBQT" \
    --center_x $CENTER_X --center_y $CENTER_Y --center_z $CENTER_Z \
    --size_x $GRID_SIZE --size_y $GRID_SIZE --size_z $GRID_SIZE \
    --exhaustiveness $EXHAUSTIVENESS \
    --out "$OUTPUT_DOCKED" \
    > "$OUTPUT_LOG" 2>&1

# --- 3. GENERATE SUMMARY ---
if [ -f "$OUTPUT_LOG" ]; then
    # Extract the top affinity score (first mode)
    # Vina log format usually has a table, we grab the first line starting with '   1 '
    AFFINITY=$(grep "   1 " "$OUTPUT_LOG" | head -n 1 | awk '{print $2}')
    
    # Get filenames for the report
    REC_NAME=$(basename "$RECEPTOR_PDBQT" .pdbqt)
    LIG_NAME=$(basename "$LIGAND_PDBQT" .pdbqt)
    
    echo "Receptor,Ligand,Affinity_kcal_mol" > "$OUTPUT_SUMMARY"
    echo "$REC_NAME,$LIG_NAME,$AFFINITY" >> "$OUTPUT_SUMMARY"
    
    echo "Docking Complete. Best Affinity: $AFFINITY"
else
    echo "Error: Vina log not found. Docking failed."
    exit 1
fi