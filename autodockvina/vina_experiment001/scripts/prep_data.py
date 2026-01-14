from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

# 1. Read the SMILES string (assuming content of your sml file)
smiles_string = open('/workspace/vina_experiment001/data/library1_smiles/Glutamine.sml').read().strip()
mol = Chem.MolFromSmiles(smiles_string)

# 2. Add Hydrogens and Generate 3D Coordinates
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)

# 3. Optional: Write to SDF if you specifically need that file
w = Chem.SDWriter('/workspace/vina_experiment001/data/library1_sdf/Glutamine.sdf')
w.write(mol)
w.close()

# 4. Prepare directly for Vina (PDBQT) using Meeko
# Meeko handles the specialized atom types Vina needs automatically
preparator = MoleculePreparation()
preparator.prepare(mol)
pdbqt_string = preparator.write_pdbqt_string()

with open('/workspace/vina_experiment001/data/library1_pdbqt/Glutamine.pdbqt', 'w') as f:
    f.write(pdbqt_string)

print("Conversion complete: ligand.sdf and ligand.pdbqt created.")