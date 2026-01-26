# python3 /workspace/workflow/scripts/bioemu_pipeline.py --fasta /workspace/data/inputs/target.fasta --working_dir data/interim/bioemu_workdir --output_pdb data/inputs/bioemu_target.pdb --samples 10
import sys
import os
import subprocess
import numpy as np
import random
import shutil
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import torch
import MDAnalysis as mda
import glob

# --- CONFIGURATION ---
AA1_TO_3 = {
    'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS',
    'I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN',
    'R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'
}

def run_bioemu(fasta_path, output_dir, num_samples=10):
    """Executes the BioEmu sampling command."""
    print(f"Starting BioEmu sampling for {fasta_path}...")
    # --- GPU DIAGNOSTIC CHECK ---
    if torch.cuda.is_available():
        device_name = torch.cuda.get_device_name(0)
        device_count = torch.cuda.device_count()
        print(f"GPU DETECTED: {device_name} (Total GPUs: {device_count})")
        print("BioEmu will run in FAST mode.")
    else:
        print("NO GPU DETECTED. BioEmu will run in SLOW mode (CPU).")
        print("   If you expected a GPU, check your Docker run command (--gpus all).")

    cmd = [
        sys.executable, "-m", "bioemu.sample",
        "--sequence", fasta_path,
        "--num_samples", str(num_samples),
        "--output_dir", output_dir
    ]
    
    # Run command and wait for it to finish
    try:
        subprocess.run(cmd, check=True)
        print("BioEmu sampling complete.")
    except subprocess.CalledProcessError as e:
        print(f"Error running BioEmu: {e}")
        exit(1)

def convert_npz_to_pdb(bioemu_workdir):
    """
    Converts BioEmu .npz samples to PDB files using the topology.pdb from the same run.
    
    Args:
        bioemu_workdir (str): Path to the BioEmu output directory 
                              (e.g., "data/interim/bioemu_workdir").
    """
    # 1. Define Paths
    topology_path = os.path.join(bioemu_workdir, "ensembles", "topology.pdb")
    npz_dir = os.path.join(bioemu_workdir, "ensembles")
    pdb_out_dir = os.path.join(bioemu_workdir, "pdb_backbone")

    print(f"Converting NPZ to PDB in {pdb_out_dir}...")
    os.makedirs(pdb_out_dir, exist_ok=True)

    # 2. Load Topology
    if not os.path.exists(topology_path):
        print(f"Error: Topology file not found at {topology_path}")
        return []

    try:
        u = mda.Universe(topology_path)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        return []

    # 3. Select Atoms (N, CA, C, O)
    # BioEmu outputs 4 backbone atoms per residue.
    # We select exactly these to match the shape of the .npz data.
    backbone = u.select_atoms("name N CA C O")
    
    print(f"Topology loaded: {len(u.residues)} residues, {backbone.n_atoms} backbone atoms.")

    # 4. Find .npz files
    npz_files = sorted(glob.glob(os.path.join(npz_dir, "batch_*.npz")))
    
    if not npz_files:
        print("No .npz files found.")
        return []

    generated_pdbs = []
    counter = 1

    # 5. Iterate and Write
    for npz_path in npz_files:
        try:
            data = np.load(npz_path, allow_pickle=True)
            
            # 'pos' shape is typically (Batch_Size, Num_Residues, 4, 3)
            # We reshape it to (Batch_Size, Num_Atoms, 3) to match MDAnalysis
            coords_batch = data["pos"]
            
            for coords in coords_batch:
                # Flatten: (Residues, 4, 3) -> (Total_Atoms, 3)
                flat_coords = coords.reshape(-1, 3)
                
                # Safety Check
                if len(flat_coords) != backbone.n_atoms:
                    print(f"Shape mismatch in sample {counter}: "
                          f"NPZ has {len(flat_coords)} atoms, Topology has {backbone.n_atoms}.")
                    continue
                
                # Update positions in the Universe
                backbone.positions = flat_coords
                
                # Write PDB
                pdb_filename = f"sample_{counter:06d}.pdb"
                out_path = os.path.join(pdb_out_dir, pdb_filename)
                
                # write_selection uses the atom names/resnames from topology.pdb
                backbone.write(out_path)
                
                generated_pdbs.append(out_path)
                counter += 1
                
        except Exception as e:
            print(f"Error processing {npz_path}: {e}")

    print(f"Successfully created {len(generated_pdbs)} PDB files in {pdb_out_dir}")
    return generated_pdbs

def pick_random_structure(pdb_list, final_output_path):
    """Picks one random PDB from the list and copies it."""
    if not pdb_list:
        print("No PDBs available to select.")
        return

    selected_pdb = random.choice(pdb_list)
    print(f"Randomly selected: {os.path.basename(selected_pdb)}")
    
    shutil.copy(selected_pdb, final_output_path)
    print(f"Saved for docking as: {final_output_path}")

def add_sidechains(input_pdb, output_pdb):
    """
    Reconstructs missing sidechains and atoms using PDBFixer.
    """
    print(f"Reconstructing sidechains for {input_pdb}...")
    
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
            
        print(f"Full atom structure saved to {output_pdb}")
        return True
        
    except Exception as e:
        print(f"Failed to reconstruct sidechains: {e}")
        return False

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BioEmu Pipeline")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA")
    parser.add_argument("--output_pdb", required=True, help="Path for the final selected PDB")
    parser.add_argument("--working_dir", required=True, help="Directory for intermediate files")
    parser.add_argument("--samples", type=int, default=10, help="Number of samples")
    
    args = parser.parse_args()

    # Define paths based on arguments
    ensembles_dir = os.path.join(args.working_dir, "ensembles")
    
    # 1. Run BioEmu
    run_bioemu(args.fasta, ensembles_dir, args.samples)

    # 2. Convert to PDBs
    all_pdbs = convert_npz_to_pdb(args.working_dir)

    if not all_pdbs:
        print("Error: No PDBs generated.")
        sys.exit(1)

    # 3. Pick ONE Random Winner
    # We save it to a temporary path first because it only has backbone atoms
    temp_backbone = os.path.join(args.working_dir, "selected_backbone_temp.pdb")
    pick_random_structure(all_pdbs, temp_backbone)

    # 4. Reconstruct Sidechains & Hydrogens
    # We take the temp backbone, fix it, and save it to the FINAL output path
    success = add_sidechains(temp_backbone, args.output_pdb)
    
    if success:
        print(f"Pipeline finished successfully. Final structure: {args.output_pdb}")
        # Optional: Clean up temp file
        if os.path.exists(temp_backbone):
            os.remove(temp_backbone)
    else:
        print("Pipeline failed at sidechain reconstruction step.")
        sys.exit(1)