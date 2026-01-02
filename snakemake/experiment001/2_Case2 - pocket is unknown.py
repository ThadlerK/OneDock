import streamlit as st
import py3Dmol
from stmol import showmol
import re
import ipywidgets
import shutil
import subprocess
import pandas as pd
from pathlib import Path, PurePath
import zipfile
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from tools.fpocket import run_fpocket


st.title('Predict your pockets from scratch')
st.write('so let\'s get started :)')
st.write('First, let\'s look at your protein structure.')
pdb_file = st.file_uploader(label = 'Upload your PBD file here',
                            type = 'pdb')

#here we visualize the pdb file
if pdb_file is not None:
    pdb_str = pdb_file.getvalue().decode("utf-8") #decodes bytes into string
    view = py3Dmol.view(width=800, height= 500) #sets the room and its size
    view.addModel(pdb_str, "pdb") #adds protein structure
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo() #zooms to protein
    showmol(view, height = 500, width = 800)

##################################################################################################
######################################### pocket detection########################################
##################################################################################################

st.header('Pocket detection')
st.write('Since you don\'t know the binding pocet of your protein, we will use two \
         different tools - fpocket anf P2Rank - to predict potential binding pockets.')

#################################parameters for pocket detection##################################

st.subheader('parameters for pocket detection')
st.write('Here, you can set the parameters for fpocket.')

col1, col2, col3 = st.columns(3)
with col1:
    fpocket_min_volume = st.select_slider(
        'minimum pocket volume',
        options=list(range(0, 1001)),
        value=400) #default position

with col2:
    fpocket_max_volume = st.select_slider(
        'maximum pocket volume',
        options=list(range(400, 3001)),
        value=1800) #default position

with col3:
    fpocket_drug_score = st.select_slider(
        'minimum druggability score',
        options=[i / 100 for i in range(0, 101)],
        value=0.5) #default position
    if fpocket_drug_score < 0.2:
        st.error('unlikely druggable pockets might be included')
    if fpocket_drug_score >= 0.2 and fpocket_drug_score < 0.4:
        st.warning('pockets with weak druggability score might be included')
    if fpocket_drug_score >= 0.4 and fpocket_drug_score < 0.7:
        st.warning('pockets are potentially druggable')
    if fpocket_drug_score >= 0.7:
        st.success('pockets are highly druggable')


col1, col2 = st.columns(2)
with col1: 
    fpocket_nr_alpha_spheres = st.select_slider(label = 'number of alpha spheres',
                                                 options = list(range(1,50)), #check for best range
                                                 value = 35)
with col2: 
    fpocket_score = st.select_slider(label = 'minimum fpocket score',
                                     options = [i / 100 for i in range(0,101)],
                                     value = 0)


st.write('And here you can set the parameters for P2Rank.')

col1, col2, col3 = st.columns(3)
with col1:
    P2rank_score = st.select_slider(label = 'minimum P2Rank score',
                                    options = list(range(0,101)),
                                    value = 60)
    if P2rank_score < 40:
        st.error('low confidence pockets might be included')
    if P2rank_score >=40 and P2rank_score <60:
        st.warning('moderate confidence pockets might be included')
    if P2rank_score >=60:
        st.success('only high confidence pockets are included')

with col2:
    P2rank_probability = st.select_slider(label = 'minimum P2Rank probability',
                                          options = [i / 100 for i in range(0,101)],
                                          value = 0.5)
with col3:
    P2rank_rank = st.select_slider(label = 'maximum P2Rank rank',
                                   options = list(range(1,51)),
                                   value = 15)
    
st.write('during pocket prediction, numerous files will be created. \
         Where should we store them for you?')
output_path = st.text_input(label = 'path for you output directory')


###########################################fpocket#################################################

st.subheader('fpocket prediction')

if st.button('run fpocket'):
    pdb_path = Path(output_path) / pdb_file.name
    with open(pdb_path, "wb") as f:
        f.write(pdb_file.getbuffer())
    st.write('running fpocket...')
    run_fpocket(
        pdb = pdb_path,
        pockets_zip = Path(output_path) / 'output_fpocket' / 'temp' / 'all_fpockets.zip',
        summary= Path(output_path) / 'output_fpocket'/ 'temp' / 'fpocket_summary.json',

        num_spheres= fpocket_nr_alpha_spheres,
        filter = True,
        filter_summary= Path(output_path) / 'output_fpocket' / 'temp' / 'filtered_fpockets.zip',
        pockets_dir = Path(output_path) / 'output_fpocket' / 'filtered_pockets',
        min_volume = fpocket_min_volume,
        max_volume = fpocket_max_volume,
        drug_score = fpocket_drug_score,
        score = fpocket_score
    )
    st.write('fpocket run complete!')


###########################################P2Rank##################################################

st.subheader('P2Rank prediction')
st.write('to run this skript, P2Rank must be installed and the path to the \
         P2Rank file must be known.')
st.write('You can download P2Rank here: https://github.com/rdk/p2rank/releases')
st.write('Please enter the directory path to the P2Rank file\
         and the directory you want the output to be stored in:')

#inputs
path_P2Rank = st.text_input(label = 'path to P2Rank file')
path_output_P2Rank = Path(output_path) / 'P2Rank_output'

if st.button("run P2Rank"):
    path_output_P2Rank.mkdir(parents=True, exist_ok=True)
    #save the pdb file
    pdb_path = Path(path_output_P2Rank) / pdb_file.name
    with open(pdb_path, "wb") as f:
        f.write(pdb_file.getbuffer())

    if path_P2Rank == '':
        st.warning('please provide path to the P2Rank file')
    if path_output_P2Rank == '':
        st.warning('please provide path to output directory for P2Rank')
    if path_P2Rank != '' and path_output_P2Rank != '':
        subprocess.run(
            [path_P2Rank, "predict", "-f", pdb_path, "-o", str(path_output_P2Rank)],
            check=True
        )

        csv_path = Path(path_output_P2Rank) / f"{pdb_file.name}_predictions.csv"

        if csv_path.exists():
            df = pd.read_csv(csv_path)
            st.success("Predictions generated")
            st.dataframe(df)
        else:
            st.error("predictions.csv not found")
        with open(csv_path, "rb") as f:
            st.download_button(
                label=f"Download {csv_path.name}",
                data=f,
                file_name=csv_path.name,
                mime="text/csv"
            )




