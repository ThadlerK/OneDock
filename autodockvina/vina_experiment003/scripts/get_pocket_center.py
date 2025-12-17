import sys
import re
import argparse

def calculate_center(pdb_file, csv_file, target_protein_name, desired_site, desired_conf):
    """
    Berechnet den Mittelpunkt basierend auf flexiblen Filtern.
    desired_site: z.B. 'orthosteric', 'allosteric' oder 'all'
    desired_conf: z.B. 'occluded', 'inward', 'outward' oder 'auto' (nimmt PDB Namen)
    """
    target_res_ids = set()
    
    # 1. CSV PARSEN UND FILTERN
    try:
        with open(csv_file, 'r') as f:
            for line in f:
                if 'ResidueSite' in line or not line.strip(): continue
                parts = line.strip().split(',')
                if len(parts) < 6: continue
                
                # CSV Spalten lesen
                residue_code = parts[0].strip()
                site_type = parts[1].strip()     # z.B. orthosteric
                csv_conf = parts[4].strip()      # z.B. occluded
                csv_prot = parts[5].strip()      # z.B. SLC6A20

                # FILTER 1: PROTEIN NAME
                if csv_prot.lower() not in target_protein_name.lower() and target_protein_name.lower() not in csv_prot.lower():
                    continue

                # FILTER 2: SITE TYPE (User Input)
                if desired_site != 'all' and site_type.lower() != desired_site.lower():
                    continue

                # FILTER 3: CONFORMATION
                # Wenn User 'auto' sagt, versuchen wir es aus dem Dateinamen zu erraten
                current_conf_filter = desired_conf
                if desired_conf == 'auto':
                    if 'occluded' in target_protein_name.lower(): current_conf_filter = 'occluded'
                    elif 'inward' in target_protein_name.lower(): current_conf_filter = 'inward'
                    elif 'outward' in target_protein_name.lower(): current_conf_filter = 'outward'
                    else: current_conf_filter = 'all' # Fallback
                
                # Check: Passt die CSV-Zeile zur gewünschten Konformation?
                # Wir nutzen "in", damit "inward-open" auf "inward" matcht
                if current_conf_filter != 'all':
                    if current_conf_filter.lower() not in csv_conf.lower():
                        continue

                # TREFFER: Residue ID extrahieren
                match = re.search(r'(\d+)', residue_code)
                if match:
                    target_res_ids.add(match.group(1))

    except Exception as e:
        sys.stderr.write(f"Error reading CSV: {e}\n")
        sys.exit(1)

    if not target_res_ids:
        # Kein Fehler, aber auch keine Koordinaten -> Leerer String an Bash zurückgeben
        # Damit weiß Bash: "Hier nichts docken."
        sys.exit(0)

    # 2. PDB KOORDINATEN HOLEN
    coords = []
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    res_seq = line[22:26].strip()
                    if res_seq in target_res_ids:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
        
        if coords:
            avg_x = sum(c[0] for c in coords) / len(coords)
            avg_y = sum(c[1] for c in coords) / len(coords)
            avg_z = sum(c[2] for c in coords) / len(coords)
            print(f"{avg_x:.3f} {avg_y:.3f} {avg_z:.3f}")
        else:
            sys.stderr.write(f"DEBUG: CSV matched residues {list(target_res_ids)[:3]}... but not found in PDB structure.\n")
            sys.exit(1)

    except Exception as e:
        sys.exit(1)

if __name__ == "__main__":
    # Wir nutzen argparse für saubere Argumente
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Path to PDB file")
    parser.add_argument("csv", help="Path to CSV file")
    parser.add_argument("name", help="Target Protein Name")
    parser.add_argument("--site", default="orthosteric", help="Desired binding site (orthosteric/allosteric/all)")
    parser.add_argument("--conf", default="auto", help="Desired conformation (occluded/inward/outward/auto)")
    
    args = parser.parse_args()
    
    calculate_center(args.pdb, args.csv, args.name, args.site, args.conf)