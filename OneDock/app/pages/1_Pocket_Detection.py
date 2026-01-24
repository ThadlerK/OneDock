import streamlit as st
import py3Dmol
from stmol import showmol
import subprocess
import pandas as pd
from pathlib import Path, PurePath
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from glob import glob
import re
import yaml 
import os
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1
from utils import filter_P2Rank, P2Rank_to_PDB, pocket_comparison, run_fpocket, save_config, load_config

# pocket_residues: 21,22,23,24,25,26,27,102,250,251,252,253,256,406,410
# ref_residues: 48,49,50,51,52,53,54,129,277,278,279,280,283,431,435

#set the session states
if "pocket_unknown" not in st.session_state:
    st.session_state.pocket_unknown = False
if "fpocket_run" not in st.session_state:
    st.session_state.fpocket_run = False
if "P2Rank_run" not in st.session_state:
    st.session_state.P2Rank_run = False
if "pocket_comparison_run" not in st.session_state:
    st.session_state.pocket_comparison_run = False


st.title("Pocket Detection")
# --- POCKET SELECTION ---
st.subheader("Binding Pocket Configuration")
pocket_status = st.radio("Is the binding pocket known?", ["Known", "Unknown"], index=0)


if pocket_status == "Known":
    st.info("Pocket is known. Please define coordinates.")
    target_residues = st.text_input("Pocket Residues of target receptor", placeholder="e.g., B:145,B:230")

    if st.button("Save Target Config"):
        # Save directly to config.yaml
        save_config({
            "pocket_residues": target_residues,
            "grid_size": 20,
            "exhaustiveness": 8
        })

    ref_path = load_config().get("ref_path")
    if ref_path:
        target_residues = st.text_input("Pocket Residues of reference receptor", placeholder="e.g., B:145,B:230")
        if st.button("Save Reference Config"):
            # Save directly to config.yaml
            save_config({
                "ref_residues": target_residues
            })


    st.session_state.pocket_unknown = False
    if st.button("Confirm & Proceed to Docking"):
        st.switch_page("pages/2_Docking.py")
    

else:
    st.warning("Pocket is unknown. You need to run prediction.")
    save_config({"pocket_known": False})
    st.session_state.pocket_unknown = True


##################################################################################################
######################################### pocket detection########################################
##################################################################################################

if st.session_state.pocket_unknown == True:
    #pdb_file = st.file_uploader(label = 'Upload your PBD file here',
     #                           type = 'pdb')
    pdb_file_path = load_config().get("receptor_path")
    with open(pdb_file_path, 'r') as f:
        pdb_file = f.read()
        
    output_path = 'data/interim/pocket_detection/'
    Path(output_path).mkdir(parents = True, exist_ok = True)

    st.write('Since you don\'t know the binding pocket of your protein, we will use two \
            different tools - fpocket and P2Rank - to predict potential binding pockets.')

    ###########################################fpocket#############################################

    with st.container():
        st.subheader('fpocket prediction')

        ##parameters 
        st.write('Here, you can set the parameters for your fpocket prediction.')

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
            fpocket_nr_alpha_spheres = st.select_slider(label = 'number of alpha spheres',
                                                        options = list(range(1,50)), #check for best range
                                                        value = 35)


        col1, col2, col3 = st.columns(3)

        with col1: 
            fpocket_score = st.select_slider(label = 'minimum fpocket score',
                                            options = [i / 100 for i in range(0,101)],
                                            value = 0)
            
        with col2:
            fpocket_drug_score = st.select_slider(
                'minimum druggability score',
                options=[i / 100 for i in range(0, 101)],
                value=0.5) #default position

        with col3:
            if fpocket_drug_score < 0.2:
                st.error('unlikely druggable pockets might be included')
            if fpocket_drug_score >= 0.2 and fpocket_drug_score < 0.4:
                st.warning('pockets with weak druggability score might be included')
            if fpocket_drug_score >= 0.4 and fpocket_drug_score < 0.7:
                st.warning('pockets are potentially druggable')
            if fpocket_drug_score >= 0.7:
                st.success('pockets are highly druggable')



        #run fpocket via skript
        if st.button('run fpocket'):
            with st.spinner("fpocket is running..."):
                run_fpocket(
                    pdb = pdb_file_path,
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
                st.session_state.fpocket_run = True #set session state

        #pocket visualization (appears when fpocket is done)
        if st.session_state.fpocket_run == True:
            col1, col2 = st.columns(2)
            with col1:
                st.write('Visualize your pockets here:')
                path_pockets = output_path + '/output_fpocket' + '/filtered_pockets'
                pocket_list = []
                for pock in glob(f"{path_pockets}/pocket*.pdb"): #creates list of all available pockets
                    match = re.search(r'pocket(\d+)\w+\.pdb', pock)
                    if match:
                        pocket_list.append("pocket" + match.group(1))

                pocket_nr = st.selectbox('select your pocket', options = pocket_list, 
                                        index = 0, placeholder = "select a pocket")
                if pocket_nr is not "select a pocket": #reads the selected pocket file
                    pocket = path_pockets + "/" + pocket_nr + "_atm.pdb"
                    with open(pocket, "r", encoding="utf-8") as f: #this is for visulization
                        pocket_file = f.read()
                    pocket = path_pockets + "/" + pocket_nr + "_atm.pdb"
                    with open(pocket, "r", encoding="utf-8") as f: #this is for text output
                        pocket_text = f.readlines()
                    if pocket_text is not None:
                        st.write(text[12:] + "\n" for text in pocket_text[5:7]) #write pocket info

            with col2:
                    #here we visualize the pdb file
                    if pdb_file is not None:
                        pdb_str = pdb_file
                        view = py3Dmol.view(width=400, height= 250) #sets the room and its size
                        view.addModel(pdb_str, "pdb") #adds protein structure
                        view.addModel(pocket_file, "pdb") #adds pocket structure
                        view.setStyle({"model": 0}, {"cartoon": {"color": "blue"}})
                        view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "red"}})
                        view.zoomTo() #zooms to protein
                        showmol(view, height = 250, width = 400)

    ###########################################P2Rank##############################################
    with st.container():
        st.subheader('P2Rank prediction')
        #inputs
        path_output_P2Rank = Path(output_path) / 'output_P2Rank'

        st.write('Here you can set the parameters for P2Rank.')

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            P2rank_probability = st.select_slider(label = 'minimum probability',
                                                options = [i / 100 for i in range(0,101)],
                                                value = 0.5)
        with col2:
            P2rank_rank = st.select_slider(label = 'maximum P2Rank rank',
                                        options = list(range(1,51)),
                                        value = 15)
        with col3:
            P2rank_score = st.select_slider(label = 'minimum P2Rank score',
                                            options = list(range(0,101)),
                                            value = 60)
        with col4:
            if P2rank_score < 40:
                st.error('low confidence pockets might be included')
            if P2rank_score >=40 and P2rank_score <60:
                st.warning('moderate confidence pockets might be included')
            if P2rank_score >=60:
                st.success('only high confidence pockets are included')

        #run P2Rank
        if st.button("run P2Rank"):
            with st.spinner("P2Rank is running"):
                path_output_P2Rank.mkdir(parents=True, exist_ok=True)

                if  path_output_P2Rank != '':
                    subprocess.run(
                        ["/opt/p2rank_2.5.1/prank", "predict", "-f", pdb_file_path, "-o", str(path_output_P2Rank)],
                        check=True
                    )
                    filter_P2Rank(
                        P2Rank_csv_path = Path(path_output_P2Rank) / f"{os.path.basename(pdb_file_path)}_predictions.csv",
                        P2Rank_score = P2rank_score,
                        P2Rank_probability = P2rank_probability,
                        P2Rank_rank = P2rank_rank,
                        P2Rank_filtered = Path(path_output_P2Rank) / f"{os.path.basename(pdb_file_path)}_filtered.csv"
                    )

                    csv_path = Path(path_output_P2Rank) / f"{os.path.basename(pdb_file_path)}_filtered.csv"
                    P2Rank_to_PDB(
                        input_csv = csv_path,
                        input_pdb = pdb_file_path, 
                        output_pdb_dir = Path(path_output_P2Rank) / 'filtered_pockets'
                    )
                    st.session_state.P2Rank_run = True

        if st.session_state.P2Rank_run == True:
            #visualization
            col1, col2 = st.columns(2)
            with col1:
                st.write('Visualize your pockets here:')
                path_pockets2 = Path(path_output_P2Rank) / 'filtered_pockets'
                pocket_list2 = []
                for pock in glob(f"{path_pockets2}/pocket*.pdb"): #creates list of all available pockets
                    match2 = re.search(r'pocket(\d+)\w+\.pdb', pock)
                    if match2:
                        pocket_list2.append("pocket" + match2.group(1))

                pocket_nr2 = st.selectbox('select your pocket', options = pocket_list2, 
                                        index = 0, placeholder = "select a pocket")
                if pocket_nr2 is not "select a pocket": #reads the selected pocket file
                    pocket2 = path_pockets2 / f"{pocket_nr2}_atm.pdb"
                    with open(pocket2, "r", encoding="utf-8") as f: #this is for visulization
                        pocket_file2 = f.read()
                    pocket2 = path_pockets2 / f"{pocket_nr2}_atm.pdb"
                    with open(pocket2, "r", encoding="utf-8") as f: #this is for text output
                        pocket_text2 = f.readlines()
                    if pocket_text2 is not None:
                        st.write(text[10:] + "\n" for text in pocket_text2[1:4]) #write pocket info

            with col2:
                    #here we visualize the pdb file
                    if pdb_file is not None:
                        pdb_str = pdb_file
                        view = py3Dmol.view(width=400, height= 250) #sets the room and its size
                        view.addModel(pdb_str, "pdb") #adds protein structure
                        view.addModel(pocket_file2, "pdb") #adds pocket structure
                        view.setStyle({"model": 0}, {"cartoon": {"color": "blue"}})
                        view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "red"}})
                        view.zoomTo() #zooms to protein
                        showmol(view, height = 250, width = 400)
    st.space("medium")

    #########################################Method comparison######################################
    with st.container():
        st.divider() #layout
        if st.session_state.fpocket_run == True and st.session_state.P2Rank_run == True:
            st.write('Now that both pocket prediction methods have been run, \
                    you might wonder whether they predict the same pockets. \
                    Therefore, we will now compare the centers of the predicted pockets.')

            col1, col2, col3 = st.columns(3)
            with col1:
                threshold = st.select_slider(label = 'minimum distance (Å) between pocket centers', 
                                            options = [i for i in range(0,21)],
                                            value = 6)
            with col2: 
                if threshold <=6:
                    st.success('only very similar pocket will be a match')
                if threshold >6 and threshold <=12:
                    st.warning('moderately similar pockets might be a match')
                if threshold >12:
                    st.error('dissimilar pockets might be a match')
            with col3:
                if st.button('compare fpocket and P2Rank'):
                    pocket_comparison(
                        fpocket_dir = f'{output_path}/output_fpocket/filtered_pockets',
                        p2Rank_csv = f'{path_output_P2Rank}/{os.path.basename(pdb_file_path)}_filtered.csv',
                        threshold = threshold,
                        output_csv = f'{output_path}pocket_comparison.csv'
                    )
                    st.session_state.pocket_comparison_run = True

        if st.session_state.pocket_comparison_run == True:
            st.write('These are the matching pockets:')
            if pd.read_csv(f'{output_path}pocket_comparison.csv').empty:
                st.warning('no matching pockets found')
            else:
                st.dataframe(pd.read_csv(f'{output_path}pocket_comparison.csv').iloc[:, 1:])

    #########################################pocket residues######################################
    
    with st.container():
        if st.session_state.fpocket_run == True or st.session_state.P2Rank_run == True:
            st.divider() #layout
            st.write('You will need the pocket residues for docking.\
                    check them here for each of your pockets')
            if st.selectbox('select your pocket detection method', 
                            options = ['fpocket', 'P2Rank'], index = 0) == 'fpocket':
                pockets_select = st.selectbox('select your fpocket', options = pocket_list, 
                                            index = 0, placeholder = "select a pocket")
                pockets_select_path = f'{path_pockets}/{pockets_select}_atm.pdb'
            else:
                pockets_select = st.selectbox('select your P2Rank pocket', options = pocket_list2, 
                                            index = 0, placeholder = "select a pocket")
                pockets_select_path = f'{path_pockets2}/{pockets_select}_atm.pdb'
            if pockets_select is not "select a pocket":
                parser = PDBParser(QUIET = True)
                structure = parser.get_structure('pocket_structure', pockets_select_path)
                residues = set()
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[0] == ' ': #excludes heteroatoms etc.
                                resname = residue.resname
                                resnum = residue.id[1]
                                chain_id = chain.id #eig nicht nötig
                                
                                try: #sicher gehen dass er nicht verkackt
                                    one_letter = protein_letters_3to1[resname.capitalize()]
                                    residues.add(f'{one_letter}{resnum}')
                                except KeyError:
                                    pass
                st.write(','.join(sorted(residues)))
                if st.button("save the residues for docking"):
                    save_config({"pocket-residues": ','.join(sorted(residues))})
                                
