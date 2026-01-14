import sys
import os
import subprocess
import numpy as np
import random
import shutil
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# --- CONFIGURATION ---
AA1_TO_3 = {
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS',
    'I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN',
    'R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'
}

def run_bioemu(fasta_path, output_dir, num_samples=100):
    """Executes the BioEmu sampling command."""
    print(f"🚀 Starting BioEmu sampling for {fasta_path}...")
    
    cmd = [
        sys.executable, "-m", "bioemu.sample",
        "--sequence", fasta_path,
        "--num_samples", str(num_samples),
        "--output_dir", output_dir
    ]
    
    # Run command and wait for it to finish
    try:
        subprocess.run(cmd, check=True)
        print("✅ BioEmu sampling complete.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running BioEmu: {e}")
        exit(1)

def convert_npz_to_pdb(npz_dir, pdb_out_dir):
    """Converts BioEmu .npz files to CA-only PDB files."""
    print(f"🔄 Converting NPZ to PDB in {pdb_out_dir}...")
    os.makedirs(pdb_out_dir, exist_ok=True)
    
    npz_files = sorted([f for f in os.listdir(npz_dir) if f.startswith("batch_") and f.endswith(".npz")])
    
    if not npz_files:
        print("⚠️ No .npz files found! Did BioEmu run correctly?")
        return []

    generated_pdbs = []
    counter = 1

    for npz_name in npz_files:
        npz_path = os.path.join(npz_dir, npz_name)
        data = np.load(npz_path, allow_pickle=True)

        coords_all = data["pos"]
        seq = data["sequence"]
        
        # Handle sequence format quirks
        if isinstance(seq, np.ndarray): seq = seq.item()
        if isinstance(seq, bytes): seq = seq.decode()

        # Iterate through every sample in the batch
        for coords in coords_all:
            pdb_filename = f"sample_{counter:06d}.pdb"
            pdb_path = os.path.join(pdb_out_dir, pdb_filename)

            with open(pdb_path, "w") as f:
                for i, (x, y, z) in enumerate(coords, start=1):
                    # Safety check for sequence length mismatch
                    if i > len(seq): break
                    
                    aa3 = AA1_TO_3.get(seq[i-1], "UNK")
                    # Standard PDB Atom Line Format
                    line = (
                        f"ATOM  {i:5d}  CA  {aa3} A{i:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
                    )
                    f.write(line + "\n")
                f.write("END\n")
            
            generated_pdbs.append(pdb_path)
            counter += 1
            
    print(f"✅ Converted {len(generated_pdbs)} PDB files.")
    return generated_pdbs

def pick_random_structure(pdb_list, final_output_path):
    """Picks one random PDB from the list and copies it."""
    if not pdb_list:
        print("❌ No PDBs available to select.")
        return

    selected_pdb = random.choice(pdb_list)
    print(f"🎲 Randomly selected: {os.path.basename(selected_pdb)}")
    
    shutil.copy(selected_pdb, final_output_path)
    print(f"💾 Saved for docking as: {final_output_path}")

def add_sidechains(input_pdb, output_pdb):
    """
    Reconstructs missing sidechains and atoms using PDBFixer.
    """
    print(f"🔧 Reconstructing sidechains for {input_pdb}...")
    
    try:
        # 1. Load the backbone structure
        fixer = PDBFixer(filename=input_pdb)
        
        # 2. Find missing residues (BioEmu might skip some)
        fixer.findMissingResidues()
        
        # 3. Find missing atoms (This detects that sidechains are missing)
        fixer.findMissingAtoms()
        
        # 4. Add the missing atoms (This builds the sidechains)
        fixer.addMissingAtoms()
        
        # 5. Add Hydrogens (Crucial for docking!)
        fixer.addMissingHydrogens(7.0) # pH 7.0
        
        # 6. Save the full structure
        with open(output_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
        print(f"✅ Full atom structure saved to {output_pdb}")
        return True
        
    except Exception as e:
        print(f"❌ Failed to reconstruct sidechains: {e}")
        return False

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BioEmu Pipeline")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA")
    
    # CHANGED: We now ask for the specific final output file path
    parser.add_argument("--output_pdb", required=True, help="Path for the final selected PDB")
    
    # We still need a directory for the intermediate files (ensembles)
    parser.add_argument("--working_dir", required=True, help="Directory for intermediate files")
    
    parser.add_argument("--samples", type=int, default=100, help="Number of samples")
    
    args = parser.parse_args()

    # Define paths based on arguments
    ensembles_dir = os.path.join(args.working_dir, "ensembles")
    pdb_dir = os.path.join(args.working_dir, "pdb_structures")

    # 1. Run BioEmu
    run_bioemu(args.fasta, ensembles_dir, args.samples)

    # 2. Convert to PDBs
    all_pdbs = convert_npz_to_pdb(ensembles_dir, pdb_dir)

    # 3. Pick Random Winner and save to the file Snakemake expects
    pick_random_structure(all_pdbs, args.output_pdb)

    temp_backbone = os.path.join(args.working_dir, "temp_backbone.pdb")

    pick_random_structure(all_pdbs, temp_backbone)

    add_sidechains(temp_backbone, args.output_pdb)