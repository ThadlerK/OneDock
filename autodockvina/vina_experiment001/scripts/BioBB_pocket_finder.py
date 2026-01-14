import sys
import os
import zipfile
import glob
import re
import logging

# --- LOGGING CONFIGURATION (CRITICAL FIX) ---
# BioBB prints logs to stdout by default. We must suppress them so they
# don't get captured by the bash script as "data".
logging.getLogger("biobb_structure_utils").setLevel(logging.ERROR)
logging.getLogger("biobb_vs").setLevel(logging.ERROR)
logging.getLogger("biobb_common").setLevel(logging.ERROR)
# Mute the root logger as well to be safe
logging.getLogger().setLevel(logging.ERROR)
# --------------------------------------------

from pathlib import Path
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

def get_pocket_center(pdb_file):
    """Calculates geometric center of a pocket PDB file."""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    # PDB Format: X (30-38), Y (38-46), Z (46-54)
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    continue
    if not coords:
        return None
    
    # Calculate geometric center (arithmetic mean of coordinates)
    avg_x = sum(c[0] for c in coords) / len(coords)
    avg_y = sum(c[1] for c in coords) / len(coords)
    avg_z = sum(c[2] for c in coords) / len(coords)
    
    return f"{avg_x:.3f} {avg_y:.3f} {avg_z:.3f}"

def run_biobb_pipeline(input_pdb, output_dir):
    # 1. Setup Paths
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    cleaned_pdb = os.path.join(output_dir, "protein_clean.pdb")
    pockets_zip = os.path.join(output_dir, "pockets.zip")
    pockets_summary = os.path.join(output_dir, "summary.json")
    filtered_zip = os.path.join(output_dir, "filtered_pockets.zip")
    
    # 2. Extract Molecule (Remove existing ligands/waters to avoid interference)
    try:
        extract_molecule(input_structure_path=input_pdb,
                         output_molecule_path=cleaned_pdb)
    except Exception as e:
        # Fallback: if extraction fails (e.g., purely protein PDB), use original
        cleaned_pdb = input_pdb

    # 3. Run fpocket (BioBB)
    # Finding cavities on the surface
    prop_run = {"min_radius": 3, "max_radius": 6, "num_spheres": 35}
    try:
        fpocket_run(input_pdb_path=cleaned_pdb,
                    output_pockets_zip=pockets_zip,
                    output_summary=pockets_summary,
                    properties=prop_run)
    except Exception as e:
        # Print to stderr so it doesn't break the bash pipe
        sys.stderr.write(f"Error running fpocket: {e}\n")
        return

    # 4. Filter Pockets (BioBB)
    # Selecting the best pockets based on druggability and volume
    prop_filter = {
        "volume": [500, 2000], 
        "druggability_score": [0.5, 1.0],
        "score": [0, 1]
    }
    
    if not os.path.exists(pockets_summary):
        sys.stderr.write("Error: fpocket failed to generate summary.\n")
        return

    try:
        fpocket_filter(input_pockets_zip=pockets_zip,
                       input_summary=pockets_summary,
                       output_filter_pockets_zip=filtered_zip,
                       properties=prop_filter)
    except Exception as e:
        sys.stderr.write(f"Error filtering pockets: {e}\n")
        return

    # 5. Extract the Filtered Pockets
    extracted_pockets_dir = os.path.join(output_dir, "pockets_extracted")
    
    if not os.path.exists(filtered_zip):
        # This usually means no pockets met the filter criteria
        return

    with zipfile.ZipFile(filtered_zip, 'r') as zip_ref:
        zip_ref.extractall(extracted_pockets_dir)

    # 6. Find PDB files and Print Coordinates
    # We use recursive globbing (rglob) in case the zip contained subfolders
    # This ensures we find the .pdb files wherever they are inside the extraction dir
    pocket_files = list(Path(extracted_pockets_dir).rglob("*.pdb"))
    
    if not pocket_files:
        sys.stderr.write("Warning: Zip extracted but no PDB files found.\n")
        return

    for p_path in pocket_files:
        filename = p_path.name
        # Extract pocket ID (e.g., 'pocket1_atm.pdb' -> 1)
        match = re.search(r'pocket(\d+)', filename)
        if match:
            p_id = match.group(1)
            center = get_pocket_center(str(p_path))
            if center:
                # Output format for Bash script: ID SCORE CENTER_X CENTER_Y CENTER_Z
                # We use a placeholder score of 0.99
                print(f"{p_id} 0.99 {center}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        # Print usage to stderr
        sys.stderr.write("Usage: python3 run_biobb_pockets.py <input_pdb> <output_dir>\n")
        sys.exit(1)
    
    run_biobb_pipeline(sys.argv[1], sys.argv[2])
    sys.stdout.flush() # Ensure buffer is flushed for Bash to read