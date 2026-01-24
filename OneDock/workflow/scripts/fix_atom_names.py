import sys

mol2_file = sys.argv[1]
pdb_in = sys.argv[2]
pdb_out = sys.argv[3]

# 1. Read Atom Names from Mol2
mol2_names = []
in_atom = False
with open(mol2_file) as f:
    for line in f:
        if line.startswith("@<TRIPOS>ATOM"): in_atom = True; continue
        if line.startswith("@<TRIPOS>"): in_atom = False
        if in_atom and line.strip():
            parts = line.split()
            if len(parts) >= 2: mol2_names.append(parts[1])

# 2. Apply to PDB
with open(pdb_in) as f_in, open(pdb_out, 'w') as f_out:
    i = 0
    for line in f_in:
        if line.startswith(("HETATM", "ATOM")) and line[17:20].strip() in ["LIG", "UNL"]:
            if i < len(mol2_names):
                # Replace atom name (cols 13-16) with name from mol2
                new_name = mol2_names[i][:4].ljust(4)
                # Force resname to LIG
                line = line[:12] + new_name + line[16:17] + "LIG" + line[20:]
                i += 1
        f_out.write(line)