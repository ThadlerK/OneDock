#This skript allows to convert filtered (or unfiltered) P2Rank csv output into pdb files 
#that can be visualized.

import pandas as pd
from pathlib import Path

def P2Rank_to_PDB(
        input_csv: str,
        input_pdb: str,
        output_pdb_dir: str
        ):
    """
    This function converts P2Rank (filtered) csv output in to pdb files, 
    that will be stored in a new directory, one file per pocket. 

    Arguments:
        input_csv: 
            path to the (filtered) P2Rank csv file
        input_pdb:
            path to the input pdb file used for P2Rank prediction
        output_pdb_dir:
            path to the output directory storing all the pdb files

    """
    #make sure the output directory exists
    Path(output_pdb_dir).mkdir(parents = True, exist_ok = True)

    #read the P2Rank csv file
    P2pockets = pd.read_csv(input_csv)

    #extract the atoms from the pdb file
    with open(input_pdb) as f:
        pdb_lines = [l for l in f if l.startswith(('HETATM', 'ATOM'))]

    #iterate over the file rows and extract the pdb files
    for i, row in P2pockets.iterrows():
        name = row['name'].strip()
        atom_ids = set(map(int, row['surf_atom_ids'].split()))

        out_path = Path(output_pdb_dir) / f"{name}_atm.pdb"

        with open(out_path, 'w') as out:
            # ---- HEADER / REMARKS ----
            #write the header information
            out.write(f"HEADER    P2RANK POCKET {name}\n")
            out.write(f"REMARK    Rank: {row['rank']}\n")
            out.write(f"REMARK    Score: {row['score']}\n")
            out.write(f"REMARK    Probability: {row['probability']}\n")

            #write the atom information
            for line in pdb_lines:
                if int(line[6:11]) in atom_ids:
                    out.write(line)        