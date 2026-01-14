#this is the fpocket script for streamlit
#it includes filtering options
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter
from pathlib import Path
import zipfile

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

