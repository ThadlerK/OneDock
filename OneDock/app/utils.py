# app/utils.py
import yaml
import os
import glob
import zipfile
import shutil
import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path
from Bio.PDB import PDBParser
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter




CONFIG_PATH = "config/config.yaml"
DATA_DIRS = [
    "data/inputs",
    "data/inputs/library_split",
    "data/interim", 
    "data/results",
]

def save_config(new_data):
    """
    Updates the config.yaml file with new keys/values 
    without overwriting existing ones.
    """
    os.makedirs("config", exist_ok=True)
    
    # Load existing config if it exists
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        config = {}

    # Update with new data
    config.update(new_data)

    # Write back to file
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(config, f)
    
    return config

def load_config():
    """loads the keays/values from the config.yaml file"""
    if not os.path.exists(CONFIG_PATH):
        return {}
    with open(CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


def reset_project():
    """
    Wipes all data and resets configuration.
    """
    # 1. Delete the entire data folder
    if os.path.exists("data"):
        try:
            shutil.rmtree("data")
        except Exception as e:
            return f"Error deleting data: {e}"

    # 2. Re-create empty structure
    for d in DATA_DIRS:
        os.makedirs(d, exist_ok=True)

    # 3. Reset Config File to empty
    with open(CONFIG_PATH, "w") as f:
        yaml.dump({}, f)

    return "Success"


#fpocket
def run_fpocket(
        pdb: str, 
        pockets_zip: str,
        summary:str,
        num_spheres: int,
        filter: bool ,
        filter_summary: str,
        pockets_dir: str,
        min_volume: int,
        max_volume: int,
        drug_score: float,
        score: float
        ):
    
    """
    Docstring for run_fpocket
    
    runs fpocket on given pdb file
    Arguments:
        pdb:
            path to input pdb file
        pocket_zip: 
            path to output pocket with file name (.zip)
            output of fpocket as zip file
        summary:
            path to output with file name (.json)
            give a summary of all pockets
        properties (num_spheres): 
            dictionary with necessary properties for fpocket
            including min and max radius and the number of alpha spheres
        filter_summary: 
            path to filter summary file (.json)
        
    """
    Path(pockets_zip).parent.mkdir(parents=True, exist_ok=True)
    Path(summary).parent.mkdir(parents=True, exist_ok=True)

    fpocket_run(input_pdb_path = pdb, 
                output_pockets_zip = pockets_zip,
                output_summary = summary,
                properties = {"min_radius":3, 
                            "max_radius": 6, 
                            "num_spheres": num_spheres})
    if filter == True:
        fpocket_filter(input_pockets_zip = pockets_zip,
               input_summary = summary,
               output_filter_pockets_zip =filter_summary,
               properties =  {"volume": [min_volume, max_volume], 
                              "druggability_score": [drug_score, 1],
                              "number_of_alpha_spheres": [num_spheres,],
                              "score": [score, 1]})
        #extract the filtered pockets and store them in a directory
        with zipfile.ZipFile(filter_summary, 'r') as zip_ref:
            zip_ref.extractall(pockets_dir)


#P2Rank filtering
def filter_P2Rank(
        P2Rank_csv_path: str,
        P2Rank_score: float,
        P2Rank_probability: float,
        P2Rank_rank: int,
        P2Rank_filtered: str
        ):
    
    """
    Docstring for filter_P2Rank

    filters csv output of P2Rank according to parameters
    Arguments:
        P2Rank_csv_path:
            path to input P2Rank csv file
        P2Rank score:
            minimum pocket score set by user
        P2Rank probability:
            minimum druggability probability set by user
        P2Rank_rank:
            minimum P2Rank rank set by user
        P2Rank_filtered: path to output file  
 
       """
    Path(P2Rank_filtered).parent.mkdir(parents = True, exist_ok = True) #set output directory
    P2pockets = pd.read_csv(P2Rank_csv_path) #read P2Rank csv file

    #strip empty spaces
    P2pockets.columns = P2pockets.columns.str.strip()
    P2pockets_filtered = P2pockets[
        (P2pockets['score'] >= P2Rank_score) &
        (P2pockets['probability'] >= P2Rank_probability) &
        (P2pockets['rank'] <= P2Rank_rank)
    ]

    P2pockets_filtered.to_csv(P2Rank_filtered)


#convert P2Rank pocket files to PDB files
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


#pocket distance comparison 
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
