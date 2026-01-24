import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# Snakemake inputs
input_smi = snakemake.input.smi
output_pdbqt = snakemake.output.pdbqt

try:
    with open(input_smi, 'r') as f:
        # Read just the SMILES string
        smiles = f.read().strip().split()[0]

    mol = Chem.MolFromSmiles(smiles)
    
    if mol:
        # --- FIX STARTS HERE ---
        # 1. Remove salts/solvents by keeping only the largest fragment
        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) > 1:
            # Sort by number of atoms and take the largest
            mol = max(fragments, key=lambda m: m.GetNumAtoms())
            print(f"Warning: Input contained multiple fragments ({smiles}). Keeping largest.")
        # --- FIX ENDS HERE ---

        mol = Chem.AddHs(mol)
        
        # generate 3D coords
        params = AllChem.ETKDG()
        params.useRandomCoords = True # Helps with small molecules sometimes
        res = AllChem.EmbedMolecule(mol, params)
        
        if res == -1:
            # Fallback if embedding fails
            AllChem.EmbedMolecule(mol, AllChem.ETKDG(), useRandomCoords=True)

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
    else:
        # Handle invalid SMILES
        print(f"Invalid SMILES: {smiles}")
        with open(output_pdbqt, 'w') as f:
            f.write("")
            
except Exception as e:
    print(f"Error converting {input_smi}: {e}")
    sys.exit(1)