import streamlit as st
import py3Dmol
from stmol import showmol
import subprocess
import pandas as pd
from pathlib import Path, PurePath
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from tools.fpocket import run_fpocket
from utils import save_config
from glob import glob
import re
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1
from tools.P2Rank_filtering import filter_P2Rank
from tools.P2Rank_to_PDB import P2Rank_to_PDB
from tools.pocket_comparison import pocket_comparison

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
    c1, c2, c3 = st.columns(3)
    cx = c1.number_input("Center X", value=0.0)
    cy = c2.number_input("Center Y", value=0.0)
    cz = c3.number_input("Center Z", value=0.0)
    
    # Save coordinates instantly
    save_config({
        "pocket_known": True,
        "center_x": cx, 
        "center_y": cy, 
        "center_z": cz
    })
    st.session_state.pocket_unknown = False
    if st.button("Confirm & Proceed to Docking ➡️"):
        st.switch_page("pages/3_Docking.py")
    

else:
    st.warning("Pocket is unknown. You need to run prediction.")
    save_config({"pocket_known": False})
    st.session_state.pocket_unknown = True


##################################################################################################
######################################### pocket detection########################################
##################################################################################################

if st.session_state.pocket_unknown == True:
    st.write('Let\'s get started with pocket detection!')
    st.write('First, let\'s look at your protein structure.')

    pdb_file = st.file_uploader(label = 'Upload your PBD file here',
                                type = 'pdb')

    #here we visualize the pdb file
    if pdb_file is not None:
        pdb_str = pdb_file.getvalue().decode("utf-8") #decodes bytes into string
        view = py3Dmol.view(width=800, height= 400) #sets the room and its size
        view.addModel(pdb_str, "pdb") #adds protein structure
        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.zoomTo() #zooms to protein
        showmol(view, height = 400, width = 800)

    st.write('while using this skript, numerous files will be created. \
            Where should we store them for you?')
    output_path = st.text_input(label = 'path for you output directory')
    
    st.header('Pocket detection')
    st.write('Since you don\'t know the binding pocket of your protein, we will use two \
            different tools - fpocket and P2Rank - to predict potential binding pockets.')

    ###########################################fpocket#############################################

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
        pdb_path = Path(output_path) / pdb_file.name #store the pdb file
        with open(pdb_path, "wb") as f:
            f.write(pdb_file.getbuffer())
        st.warning('running fpocket...')
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
        st.success('fpocket run complete!')
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
                    pdb_str = pdb_file.getvalue().decode("utf-8") #decodes bytes into string
                    view = py3Dmol.view(width=400, height= 250) #sets the room and its size
                    view.addModel(pdb_str, "pdb") #adds protein structure
                    view.addModel(pocket_file, "pdb") #adds pocket structure
                    view.setStyle({"model": 0}, {"cartoon": {"color": "blue"}})
                    view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "red"}})
                    view.zoomTo() #zooms to protein
                    showmol(view, height = 250, width = 400)

    ###########################################P2Rank##############################################

    st.subheader('P2Rank prediction')
    st.write('to run this skript, P2Rank must be installed and the path to the \
            P2Rank file must be known.')
    st.write('You can download P2Rank here: https://github.com/rdk/p2rank/releases')
    st.write('Please enter the directory path to the P2Rank file\
            and the directory you want the output to be stored in:')

    #inputs
    path_output_P2Rank = Path(output_path) / 'P2Rank_output'

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
                                        value = 10)
    with col4:
        if P2rank_score < 40:
            st.error('low confidence pockets might be included')
        if P2rank_score >=40 and P2rank_score <60:
            st.warning('moderate confidence pockets might be included')
        if P2rank_score >=60:
            st.success('only high confidence pockets are included')

    #run P2Rank
    if st.button("run P2Rank"):
        path_output_P2Rank.mkdir(parents=True, exist_ok=True)
        #save the pdb file
        pdb_path = Path(path_output_P2Rank) / pdb_file.name
        with open(pdb_path, "wb") as f:
            f.write(pdb_file.getbuffer())

        if  path_output_P2Rank != '':
            subprocess.run(
                ["/opt/p2rank_2.5.1/prank", "predict", "-f", pdb_path, "-o", str(path_output_P2Rank)],
                check=True
            )
            filter_P2Rank(
                P2Rank_csv_path = Path(path_output_P2Rank) / f"{pdb_file.name}_predictions.csv",
                P2Rank_score = P2rank_score,
                P2Rank_probability = P2rank_probability,
                P2Rank_rank = P2rank_rank,
                P2Rank_filtered = Path(path_output_P2Rank) / f"{pdb_file.name}_filtered.csv"
            )

            csv_path = Path(path_output_P2Rank) / f"{pdb_file.name}_filtered.csv"
            P2Rank_to_PDB(
                input_csv = csv_path,
                input_pdb = pdb_path, 
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
                    pdb_str = pdb_file.getvalue().decode("utf-8") #decodes bytes into string
                    view = py3Dmol.view(width=400, height= 250) #sets the room and its size
                    view.addModel(pdb_str, "pdb") #adds protein structure
                    view.addModel(pocket_file2, "pdb") #adds pocket structure
                    view.setStyle({"model": 0}, {"cartoon": {"color": "blue"}})
                    view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "red"}})
                    view.zoomTo() #zooms to protein
                    showmol(view, height = 250, width = 400)


    #########################################Method comparison######################################
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
                    p2Rank_csv = f'{path_output_P2Rank}/{pdb_file.name}_filtered.csv',
                    threshold = threshold,
                    output_csv = f'{output_path}pocket_comparison.csv'
                )
                st.session_state.pocket_comparison_run = True

    if st.session_state.pocket_comparison_run == True:
        st.write('These are the matching pockets:')
        if pd.read_csv(f'{output_path}pocket_comparison.csv').empty:
            st.write('no matching pockets found')
        else:
            st.dataframe(pd.read_csv(f'{output_path}pocket_comparison.csv').iloc[:, 1:])
    
    #########################################pocket residues######################################
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
                            
    #we want this formal A:145,A:230,S:450