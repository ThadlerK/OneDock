# app/utils.py
import glob
import os
import re
import requests
import shutil
import time
import yaml
import zipfile
import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path
from Bio.PDB import PDBParser
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter
from bs4 import BeautifulSoup
from typing import List, Dict, Optional

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

#ADME tools
"""
SwissADME Client for automated ADME property prediction
This module provides functions to interact with SwissADME web service
"""

class SwissADMEClient:
    """
    Client for interacting with SwissADME web service
    
    Note: SwissADME doesn't have an official API, so this uses web scraping.
    Use responsibly and respect their terms of service.
    """
    
    BASE_URL = "http://www.swissadme.ch"
    SUBMIT_URL = f"{BASE_URL}/index.php"
    
    def __init__(self, delay: float = 2.0):
        """
        Initialize the client
        
        Args:
            delay: Delay between requests in seconds (be respectful!)
        """
        self.delay = delay
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36'
        })
    
    def submit_smiles(self, smiles_list: List[str], names: Optional[List[str]] = None) -> Optional[str]:
        """
        Submit SMILES to SwissADME
        
        Args:
            smiles_list: List of SMILES strings
            names: Optional list of compound names (must match length of smiles_list)
        
        Returns:
            URL or identifier for results retrieval, or None if failed
        """
        if names is None:
            names = [f"Compound_{i+1}" for i in range(len(smiles_list))]
        
        if len(smiles_list) != len(names):
            raise ValueError("Length of smiles_list and names must match")
        
        # Format input: "SMILES Name" per line
        smiles_input = "\n".join([f"{smi} {name}" for smi, name in zip(smiles_list, names)])
        
        try:
            # Submit form
            data = {
                'smiles': smiles_input,
            }
            
            response = self.session.post(self.SUBMIT_URL, data=data, timeout=30)
            
            if response.status_code == 200:
                # Extract result URL or job ID from response
                # This depends on SwissADME's actual response format
                return response.url
            else:
                print(f"Submission failed with status code: {response.status_code}")
                return None
                
        except Exception as e:
            print(f"Error during submission: {e}")
            return None
    
    def parse_results(self, html_content: str) -> Optional[pd.DataFrame]:
        """
        Parse SwissADME HTML results into a DataFrame
        
        Args:
            html_content: HTML content from SwissADME results page
        
        Returns:
            DataFrame with ADME properties, or None if parsing failed
        """
        try:
            soup = BeautifulSoup(html_content, 'html.parser')
            
            # Find the results table (this is a placeholder - actual parsing depends on site structure)
            table = soup.find('table', {'class': 'results'})
            
            if table:
                df = pd.read_html(str(table))[0]
                return df
            else:
                print("Could not find results table in HTML")
                return None
                
        except Exception as e:
            print(f"Error parsing results: {e}")
            return None


def calculate_lipinski_violations(mol_weight: float, logp: float, h_donors: int, h_acceptors: int) -> int:
    """
    Calculate Lipinski's Rule of Five violations
    
    Args:
        mol_weight: Molecular weight (Da)
        logp: Octanol-water partition coefficient
        h_donors: Number of hydrogen bond donors
        h_acceptors: Number of hydrogen bond acceptors
    
    Returns:
        Number of violations (0-4)
    """
    violations = 0
    
    if mol_weight > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if h_donors > 5:
        violations += 1
    if h_acceptors > 10:
        violations += 1
    
    return violations


def calculate_qed(properties: Dict) -> float:
    """
    Calculate Quantitative Estimate of Drug-likeness (QED)
    Simplified version based on common molecular descriptors
    
    Args:
        properties: Dictionary with molecular properties
    
    Returns:
        QED score (0-1, higher is better)
    """
    # This is a simplified placeholder
    # Real QED calculation is more complex
    score = 0.5
    
    if 'MW' in properties and 150 <= properties['MW'] <= 500:
        score += 0.1
    if 'LogP' in properties and -0.4 <= properties['LogP'] <= 5.6:
        score += 0.1
    if 'HBA' in properties and properties['HBA'] <= 10:
        score += 0.1
    if 'HBD' in properties and properties['HBD'] <= 5:
        score += 0.1
    if 'TPSA' in properties and properties['TPSA'] <= 140:
        score += 0.1
    
    return min(1.0, score)


def format_smiles_for_swissadme(smiles_df: pd.DataFrame, smiles_col: str = 'SMILES', name_col: str = 'Name') -> str:
    """
    Format a DataFrame with SMILES for SwissADME input
    
    Args:
        smiles_df: DataFrame containing SMILES and names
        smiles_col: Name of column containing SMILES strings
        name_col: Name of column containing compound names
    
    Returns:
        Formatted string ready for SwissADME input
    """
    lines = []
    for idx, row in smiles_df.iterrows():
        smiles = row[smiles_col]
        name = row[name_col] if name_col in row else f"Compound_{idx}"
        lines.append(f"{smiles} {name}")
    
    return "\n".join(lines)


def filter_drug_like_compounds(df: pd.DataFrame, rules: str = 'lipinski') -> pd.DataFrame:
    """
    Filter compounds based on drug-likeness rules
    
    Args:
        df: DataFrame with ADME properties
        rules: Which rule set to apply ('lipinski', 'veber', 'ghose')
    
    Returns:
        Filtered DataFrame
    """
    filtered = df.copy()
    
    if rules.lower() == 'lipinski':
        if 'Lipinski #violations' in filtered.columns:
            filtered = filtered[filtered['Lipinski #violations'] == 0]
        elif all(col in filtered.columns for col in ['MW', 'LogP', 'HBD', 'HBA']):
            filtered = filtered[
                (filtered['MW'] <= 500) &
                (filtered['LogP'] <= 5) &
                (filtered['HBD'] <= 5) &
                (filtered['HBA'] <= 10)
            ]
    
    elif rules.lower() == 'veber':
        # Veber rules: TPSA ≤ 140 Ų and ≤ 10 rotatable bonds
        if all(col in filtered.columns for col in ['TPSA', 'Rotatable Bonds']):
            filtered = filtered[
                (filtered['TPSA'] <= 140) &
                (filtered['Rotatable Bonds'] <= 10)
            ]
    
    elif rules.lower() == 'ghose':
        # Ghose rules
        if all(col in filtered.columns for col in ['MW', 'LogP', 'MR', 'Atoms']):
            filtered = filtered[
                (filtered['MW'] >= 160) & (filtered['MW'] <= 480) &
                (filtered['LogP'] >= -0.4) & (filtered['LogP'] <= 5.6) &
                (filtered['MR'] >= 40) & (filtered['MR'] <= 130) &
                (filtered['Atoms'] >= 20) & (filtered['Atoms'] <= 70)
            ]
    
    return filtered


def analyze_adme_profile(df: pd.DataFrame, ligand_id: str) -> Dict:
    """
    Generate comprehensive ADME profile for a single compound
    
    Args:
        df: DataFrame with ADME properties
        ligand_id: ID of the ligand to analyze
    
    Returns:
        Dictionary with ADME profile and recommendations
    """
    compound = df[df['Name'] == ligand_id].iloc[0] if 'Name' in df.columns else None
    
    if compound is None:
        return {'error': 'Compound not found'}
    
    profile = {
        'ligand_id': ligand_id,
        'drug_likeness': {},
        'pharmacokinetics': {},
        'warnings': [],
        'score': 0
    }
    
    # Drug-likeness assessment
    if 'Lipinski #violations' in compound:
        violations = compound['Lipinski #violations']
        profile['drug_likeness']['lipinski_violations'] = int(violations)
        profile['score'] += (1 - violations/4) * 25
        
        if violations > 1:
            profile['warnings'].append(f"Multiple Lipinski violations ({violations})")
    
    # Bioavailability
    if 'Bioavailability Score' in compound:
        bioavail = compound['Bioavailability Score']
        profile['pharmacokinetics']['bioavailability_score'] = float(bioavail)
        profile['score'] += bioavail * 25
        
        if bioavail < 0.55:
            profile['warnings'].append("Low oral bioavailability predicted")
    
    # GI absorption
    if 'GI absorption' in compound:
        gi = compound['GI absorption']
        profile['pharmacokinetics']['gi_absorption'] = gi
        if gi == 'High':
            profile['score'] += 25
        
        if gi == 'Low':
            profile['warnings'].append("Low GI absorption predicted")
    
    # BBB permeability
    if 'BBB permeant' in compound:
        bbb = compound['BBB permeant']
        profile['pharmacokinetics']['bbb_permeant'] = bbb
    
    # PAINS alerts
    if 'PAINS' in compound:
        pains = compound['PAINS']
        if pains > 0:
            profile['warnings'].append(f"PAINS alert: {pains} problematic substructures detected")
            profile['score'] -= 10
    
    # Synthetic accessibility
    if 'Synthetic accessibility' in compound:
        sa_score = compound['Synthetic accessibility']
        profile['drug_likeness']['synthetic_accessibility'] = float(sa_score)
        
        if sa_score > 6:
            profile['warnings'].append("Difficult to synthesize (SA score > 6)")
    
    return profile


if __name__ == "__main__":
    # Example usage
    print("SwissADME Client Module")
    print("This module provides utilities for ADME property prediction")
    
    # Example: Calculate Lipinski violations
    example_props = {
        'mol_weight': 450,
        'logp': 3.5,
        'h_donors': 2,
        'h_acceptors': 6
    }
    
    violations = calculate_lipinski_violations(
        example_props['mol_weight'],
        example_props['logp'],
        example_props['h_donors'],
        example_props['h_acceptors']
    )
    
    print(f"\nExample: Lipinski violations = {violations}")

