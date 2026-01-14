import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# Snakemake inputs
input_smi = snakemake.input.smi
output_pdbqt = snakemake.output.pdbqt

try:
    with open(input_smi, 'r') as f:
        smiles = f.read().strip().split()[0]

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]
        
        with open(output_pdbqt, 'w') as f:
            f.write(pdbqt_string)
    else:
        # Create empty file to prevent Snakemake error, or handle gracefully
        with open(output_pdbqt, 'w') as f:
            f.write("")
except Exception as e:
    print(f"Error converting {input_smi}: {e}")
    sys.exit(1)