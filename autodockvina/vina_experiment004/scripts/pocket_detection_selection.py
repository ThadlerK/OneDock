import sys
import os
import shutil
import zipfile
import pandas as pd
import numpy as np
import glob
from pathlib import Path
from Bio.PDB import PDBParser

# Hinweis: biobb Bibliotheken müssen installiert sein. 
# Falls nicht, können Teile wie 'extract_molecule' durch einfache Biopython-Logik ersetzt werden.
try:
    from biobb_structure_utils.utils.extract_molecule import extract_molecule
    from biobb_vs.fpocket.fpocket_run import fpocket_run
    from biobb_vs.fpocket.fpocket_filter import fpocket_filter
except ImportError:
    print("BIOBB Libraries not found. Please install biobb_structure_utils and biobb_vs or run in correct environment.")
    # Fallback logik wäre nötig, aber wir gehen davon aus, dass die Umgebung stimmt.

# --- CONFIGURATION ---
BASE_DIR = Path("/workspace/vina_experiment004")
STRUCTURES_DIR = BASE_DIR / "data/structures"
P2RANK_OUTPUT_DIR = BASE_DIR / "data/p2rank_output" # Wo die p2rank CSVs liegen
OUTPUT_DIR = BASE_DIR / "data/pockets"

# Input Files (Beispiel für EINEN Rezeptor - man könnte dies auch als Loop gestalten)
# Wir nehmen den Rezeptor aus deinem Bash-Skript
RECEPTOR_NAME = "SLC6A20_8I91_complete_occluded"
PDB_FILE = STRUCTURES_DIR / f"{RECEPTOR_NAME}.pdb"
P2RANK_FILE = P2RANK_OUTPUT_DIR / "SIT1_Model_OO_no_ligands.pdb_predictions.csv"

# Settings
LIGAND_PRESENT = False # Set to True to strip HETATM
THRESHOLD_DISTANCE = 6.0 # Max distance between centers (Angstrom)

# Parameters for filtering
FPOCKET_PARAMS = {
    "volume": [500, 2500],  # Groß genug für Drug-Like
    "druggability_score": [0.3, 1.0], 
    "score": [0, 1]
}
PROPERTIES_FPOCKET = {"min_radius": 3, "max_radius": 6, "num_spheres": 30}

P2_SCORE_MIN = 5
P2_PROB_MIN = 0.5
P2_RANK_MAX = 10

# --- SETUP DIRECTORIES ---
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR = OUTPUT_DIR / "temp"
TEMP_DIR.mkdir(parents=True, exist_ok=True)
POCKETS_EXTRACT_DIR = TEMP_DIR / "pockets_extracted"
POCKETS_EXTRACT_DIR.mkdir(parents=True, exist_ok=True)

# Output Files
OUTPUT_FPOCKET_ZIP = TEMP_DIR / "all_fpockets.zip"
OUTPUT_FPOCKET_SUMMARY = TEMP_DIR / "fpocket_summary.json"
OUTPUT_FILTERED_ZIP = TEMP_DIR / "fpocket_filtered.zip"
OUTPUT_FINAL_CSV = OUTPUT_DIR / "pockets_final.csv"

# =============================================================================
# 1. PRE-PROCESSING (Ligand Removal)
# =============================================================================
print(f"--- Processing {RECEPTOR_NAME} ---")

# Wir arbeiten auf einer Kopie im Temp-Ordner, um das Original nicht zu ändern
working_pdb = TEMP_DIR / f"{RECEPTOR_NAME}_clean.pdb"

if LIGAND_PRESENT:
    print("removing ligands/heteroatoms...")
    # Alternative extraction logic if biobb fails or simple cleanup is preferred:
    # Uses Biopython to save only ATOM records (Protein)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', PDB_FILE)
    
    # Save standard PDB (Standard Bio.PDB.PDBIO handles Select automatically mostly)
    from Bio.PDB import PDBIO, Select
    class NotHetatm(Select):
        def accept_residue(self, residue):
            # 0 means HETATM/Ligand/Water in PDB ID usually
            if residue.id[0] != " ": 
                return 0
            else:
                return 1
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(working_pdb), select=NotHetatm())
else:
    shutil.copy(PDB_FILE, working_pdb)

# =============================================================================
# 2. RUN FPOCKET & FILTER
# =============================================================================
print("Running fpocket...")
fpocket_run(input_pdb_path=str(working_pdb), 
            output_pockets_zip=str(OUTPUT_FPOCKET_ZIP),
            output_summary=str(OUTPUT_FPOCKET_SUMMARY),
            properties=PROPERTIES_FPOCKET)

print("Filtering fpockets...")
fpocket_filter(input_pockets_zip=str(OUTPUT_FPOCKET_ZIP),
               input_summary=str(OUTPUT_FPOCKET_SUMMARY),
               output_filter_pockets_zip=str(OUTPUT_FILTERED_ZIP),
               properties=FPOCKET_PARAMS)

# Extract
with zipfile.ZipFile(str(OUTPUT_FILTERED_ZIP), 'r') as zip_ref:
    zip_ref.extractall(str(POCKETS_EXTRACT_DIR))

# =============================================================================
# 3. ANALYZE FPOCKETS (Calculate Centers)
# =============================================================================
fpocket_files = glob.glob(str(POCKETS_EXTRACT_DIR / "pocket*_atm.pdb"))
fpocket_data = []

print(f"Found {len(fpocket_files)} fpockets after filtering.")

for fp_file in fpocket_files:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pocket', fp_file)
    coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
    
    if len(coords) == 0: continue
    
    center = coords.mean(axis=0)
    
    # ID cleaning (pocket1_atm.pdb -> 1)
    fname = os.path.basename(fp_file)
    # Extract number using regex or simple split
    try:
        pocket_id = fname.split('_')[0].replace('pocket', '')
    except:
        pocket_id = fname

    fpocket_data.append({
        'pocket_id': pocket_id,
        'center_x': center[0],
        'center_y': center[1],
        'center_z': center[2],
        'volume': len(coords) # Approximation
    })

fpocket_df = pd.DataFrame(fpocket_data)

# =============================================================================
# 4. LOAD P2RANK PREDICTIONS
# =============================================================================
if not P2RANK_FILE.exists():
    print(f"Error: P2Rank file not found at {P2RANK_FILE}")
    print("Please run p2rank prediction first.")
    sys.exit(1)

print("Loading p2rank predictions...")
p2_df = pd.read_csv(P2RANK_FILE, skipinitialspace=True)
# Columns often have spaces: " center_x" -> "center_x"
p2_df.columns = p2_df.columns.str.strip()

# Filter P2Rank
p2_filtered = p2_df[
    (p2_df['score'] >= P2_SCORE_MIN) &
    (p2_df['probability'] >= P2_PROB_MIN) &
    (p2_df['rank'] <= P2_RANK_MAX)
].copy()

print(f"P2Rank pockets after filtering: {len(p2_filtered)}")

# =============================================================================
# 5. CONSENSUS (Distance Matching)
# =============================================================================
def calc_dist(p1, p2):
    return np.sqrt(np.sum((np.array(p1) - np.array(p2))**2))

matches = []

for idx_p2, row_p2 in p2_filtered.iterrows():
    p2_center = [row_p2['center_x'], row_p2['center_y'], row_p2['center_z']]
    
    for idx_fp, row_fp in fpocket_df.iterrows():
        fp_center = [row_fp['center_x'], row_fp['center_y'], row_fp['center_z']]
        
        dist = calc_dist(p2_center, fp_center)
        
        if dist <= THRESHOLD_DISTANCE:
            matches.append({
                'Receptor': RECEPTOR_NAME,  # <-- WICHTIG für dein Bash Script
                'p2rank_rank': row_p2['rank'],
                'fpocket_nr': row_fp['pocket_id'],
                'distance': round(dist, 2),
                # Speichere Center als formatierten String für das Bash Script parsing
                'fpocket_center': f"{fp_center[0]:.3f} {fp_center[1]:.3f} {fp_center[2]:.3f}",
                'p2Rank_center': f"{p2_center[0]:.3f} {p2_center[1]:.3f} {p2_center[2]:.3f}"
            })

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================
matches_df = pd.DataFrame(matches)

# Append to file if exists (for multi-receptor workflows), or create new
if OUTPUT_FINAL_CSV.exists():
    # Load existing to avoid duplicates if re-running
    existing_df = pd.read_csv(OUTPUT_FINAL_CSV)
    # Remove old entries for this receptor to cleanly update
    existing_df = existing_df[existing_df['Receptor'] != RECEPTOR_NAME]
    final_df = pd.concat([existing_df, matches_df], ignore_index=True)
else:
    final_df = matches_df

final_df.to_csv(OUTPUT_FINAL_CSV, index=False)

print(f"Done. Found {len(matches)} consensus pockets.")
print(f"Results saved to: {OUTPUT_FINAL_CSV}")


### Kurzanleitung zur Verwendung

# 1.  Stelle sicher, dass `p2rank` für deinen Rezeptor bereits gelaufen ist und die CSV im Ordner `data/p2rank_output` liegt.
# 2.  Lass dieses Python-Script laufen. Es erstellt (oder aktualisiert) `data/pockets/pockets_final.csv`.
# 3.  Starte danach dein `redock_specificity_v3.sh`. Es wird die neue CSV lesen, die Spalte `Receptor` finden und das Cross-Docking starten.

# ### Erklärung zur Reinigung (Ligand/ACE)
# Im Skript nutze ich `Biopython` mit der Klasse `Select(NotHetatm)`. Das ist robuster als `extract_molecule` von biobb, falls die Bibliotheken fehlen, und es entfernt zuverlässig:
# * Liganden (HETATM)
# * Wasser (HETATM)
# * Ionen

# Wenn dein ACE ein separates Protein (Kette B, C...) ist und du es auch entfernen willst, müsstest du die `accept_residue` Klasse erweitern, z.B.:
# python
# if residue.get_parent().id == "A": # Behalte nur Chain A
#     return 1
# else:
#     return 0