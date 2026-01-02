import streamlit as st

st.title("Protein Viewer")
st.write("Hier kommt die Visualisierung")
import py3Dmol
from stmol import showmol
from pathlib import Path
from biobb_vs.fpocket.fpocket_run import fpocket_run

pdb_file = st.file_uploader(label = 'Upload your PBD file here',
                            type = 'pdb')
pocket = st.file_uploader(label = 'Upload your pocket pdb file here',
                            type = 'pdb')
pocket_str = pocket.getvalue().decode("utf-8")

path_pockets = '/home/sylviane/meet-eu/Heidelberg_Team_2-1/test/output_fpocket/filtered_pockets'

pocket_list = []
for pock in path_pockets:
    g = re.findall('(?:pocket)(\d+)(?:\w+)\.(\w+)', pock)
    pocket_list.append(g[0][0])

st.selectbox('select the pocket you want to visualize', options = pocket_list, index = None)

#here we visualize the pdb file
if pdb_file is not None:
    pdb_str = pdb_file.getvalue().decode("utf-8") #decodes bytes into string
    view = py3Dmol.view(width=800, height= 500) #sets the room and its size
    view.addModel(pdb_str, "pdb") #adds protein structure
    view.addModel(pocket_str, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
    view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "red"}})
    view.zoomTo() #zooms to protein
    showmol(view, height = 500, width = 800)



