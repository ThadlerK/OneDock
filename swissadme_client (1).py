"""
SwissADME Client for automated ADME property prediction
This module provides functions to interact with SwissADME web service
"""

import requests
import time
from typing import List, Dict, Optional
import pandas as pd
from bs4 import BeautifulSoup
import re


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
