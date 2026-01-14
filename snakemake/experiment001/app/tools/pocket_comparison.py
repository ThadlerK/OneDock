#this skript is to compare fpocket and P2Rank predictions by computing the distance between the pockets
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np
import glob

def pocket_comparison(
    fpocket_dir = str,
    p2Rank_csv = str,
    threshold = float,
    output_csv = str
    ):

    """
    Docksting of pocket_comparison:
    This function compares the pockets predicted by fpocket and P2Rank by computing the distance 
    between the pocket centers. If the distance is below a certain threshold, the pockets are
    considered matching
    
    Arguments:
        fpocket_dir:
            path to folder containing fpocket pdb files
        p2Rank_csv: 
            path to the (filtered) P2Rank csv file
        threshold:
            maximum distance between pocket centers to be considered a match
        output_csv:
            path to the output csv file storing the matches

    """


    #the fpocket files first need to be converted into a common csv file:
    pdb_files = glob.glob(f"{fpocket_dir}/pocket*_atm.pdb")

    fpocket_list = []

    for pdb_file in pdb_files:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pocket', pdb_file)
        coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
        center = coords.mean(axis=0)  # center of mass of pocket atoms
        volume = len(coords)  # rough proxy for size
        pocket_id = pdb_file.split('/')[-1].replace('_atm.pdb','')
        fpocket_list.append({
            'pocket_id': pocket_id,
            'center_x': center[0],
            'center_y': center[1],
            'center_z': center[2],
            'volume': volume
        })

    fpocket_df = pd.DataFrame(fpocket_list)
    fpocket_df.to_csv(f"{fpocket_dir}/fpocket_centers.csv", index=False)


    #needed files:
    p2rank = pd.read_csv(p2Rank_csv)
    fpocket = pd.read_csv(f"{fpocket_dir}/fpocket_centers.csv")

    #we are comparing pockets via the distance of their centers. The centers need to be closer 
    #than 6 armstrong to be considered one
    #Here's a function that computes the distance between two points:
    #(d=√((x_2-x_1)²+(y_2-y_1)²+(z_2-z_1)²)) with p = (x, y, z):

    def dist(p1, p2):
        return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

    #now, let's look for matches:
    matches = []
    for i, p2 in p2rank.iterrows():
        p2_center = [p2['center_x'], p2['center_y'], p2['center_z']]
        for j, fp in fpocket.iterrows():
            fp_center = [fp['center_x'], fp['center_y'], p2['center_z']]
            if dist(p2_center, fp_center) <= threshold:
                matches.append({
                    'p2rank_rank': p2['rank'],
                    'fpocket_nr': fp['pocket_id'],
                    'distance': dist(p2_center, fp_center)
                })
    matches_df = pd.DataFrame(matches)
    matches_df.to_csv(output_csv)

