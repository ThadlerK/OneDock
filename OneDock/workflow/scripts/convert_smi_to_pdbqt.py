import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# Snakemake inputs
input_smi = snakemake.input.smi
output_pdbqt = snakemake.output.pdbqt

try:
    with open(input_smi, 'r') as f:
        # Read just the SMILES string (first token)
        smiles = f.read().strip().split()[0]

    mol = Chem.MolFromSmiles(smiles)
    
    if mol:
        # --- 1. UNSUPPORTED ELEMENT CHECK ---
        # Standard organic set + halogens: H, C, N, O, F, P, S, Cl, Br, I
        allowed_atoms = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53} 
        has_bad_atom = False

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in allowed_atoms:
                print(f"Skipping {smiles}: Contains unsupported element {atom.GetSymbol()}")
                has_bad_atom = True
                break

        if has_bad_atom:
            # Write empty file to satisfy Snakemake and exit cleanly
            with open(output_pdbqt, 'w') as f: f.write("")
            sys.exit(0) 

        # --- 2. DESALTING ---
        # Remove salts/solvents by keeping only the largest fragment
        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) > 1:
            mol = max(fragments, key=lambda m: m.GetNumAtoms())
            print(f"Warning: Input contained multiple fragments ({smiles}). Keeping largest.")

        # --- 3. 3D GENERATION ---
        mol = Chem.AddHs(mol)
        
        # Attempt 1: ETKDG with Random Coords (usually works)
        params = AllChem.ETKDG()
        params.useRandomCoords = True
        res = AllChem.EmbedMolecule(mol, params)
        
        # Attempt 2: Fallback (Try again if failed)
        if res == -1:
            print(f"Warning: First embedding failed for {smiles}, retrying...")
            params = AllChem.ETKDG()
            params.useRandomCoords = True  # Set the option inside the object
            res = AllChem.EmbedMolecule(mol, params)

        # --- CRITICAL FIX: Check if we actually got coordinates ---
        if res == -1 or mol.GetNumConformers() == 0:
            print(f"Error: Could not generate 3D coordinates for {smiles} after retries.")
            with open(output_pdbqt, 'w') as f: f.write("")
            sys.exit(0) # Exit success so pipeline continues, but file is empty

        # --- 4. MEEKO PREPARATION ---
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
            
    else:
        # Handle invalid SMILES (parsing failed)
        print(f"Invalid SMILES: {smiles}")
        with open(output_pdbqt, 'w') as f:
            f.write("")
            
except Exception as e:
    # Catch unexpected crashes (e.g. Meeko errors)
    print(f"Error converting {input_smi}: {e}")
    sys.exit(1) # Fail the job so you notice something is wrong