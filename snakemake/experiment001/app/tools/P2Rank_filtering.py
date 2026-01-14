#This is the script to filter P2Rank pockets
#before running this, P2Rank must be applied to a .pdb file like this:
#/home/sylviane/meet-eu/Heidelberg_Team_2-1/P2Rank/p2rank_2.5.1/prank predict -f /home/sylviane/meet-eu/Heidelberg_Team_2-1/fpocket_data/struc
#tures/8I92_ACE2.pdb

#needed packages
import pandas as pd
import numpy as np
from pathlib import Path

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