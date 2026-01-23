# app/Home.py
#streamlit run app/Home.py --server.address=0.0.0.0 --server.port=8501

import streamlit as st
import os
from utils import save_config
import py3Dmol
from stmol import showmol
import subprocess
from utils import reset_project
import time

st.set_page_config(
    page_title = "OneDock",
    page_icon = "app/images/logo_protein.png"
)

st.image("app/images/logo_dunkel.png")
st.markdown("<p style='text-align: center; font-weight: bold; font-size: 16px;'>The easiest way to find the most promising Ligands for your Protein</p>", 
            unsafe_allow_html = True)
st.space("small")

st.write("Welcome to your protein-ligand matchmaker! \n" \
         "If you want to find out which ligands might bind to your protein, \n" \
        "you have come to the right place. \n" \
        "OneDock can be applied in multiple ways:")


st.space("medium")
col1, col2, col3, col4, col5, col6, col7 = st.columns ([1,3,2,3,1,3,1])
with col2:
    st.image("app/images/specificity.svg")
with col4:
    st.image("app/images/unknown_sites.svg")
with col6:
    st.space("medium")
    st.image("app/images/protein_sequence.svg")
    st.space("medium")

col1, col2, col3 = st.columns(3)
with col1:
    st.write("<p style = 'color: #AB7FE1; text-align: center;'>Find the Ligands with the highest specificity for your Protein</p>", 
             unsafe_allow_html = True)

with col2:
    st.write("<p style = 'color: #AB7FE1; text-align: center;'>Identify unknown binding sites</p>",
             unsafe_allow_html = True)


with col3:
    st.write("<p style = 'color: #AB7FE1; text-align: center;'>Work with Proteins without a defined Structure</p>",
             unsafe_allow_html = True)

st.space('large')

st.write('To do this you will be guided through the following pipeline:')

st.image('app/images/pipeline.png')

st.space('large')

st.write('Once you have completed the pipeline OneDock will provide you with \n' \
         'extensive Details about your ligands:')

col1, col2, col3, col4, col5, col6, col7 = st.columns([0.5, 3, 1.5, 3, 1.5, 3, 1])

#example visualization of pockets
with col2:
    with open("app/images/example_protein.pdb", "r") as f:
        exemp_protein = f.read()
    with open("app/images/example_pocket.pdb", "r") as f:
        exemp_pocket = f.read()
    view = py3Dmol.view(width=200, height= 200) #sets the room and its size
    view.addModel(exemp_protein, "pdb") #adds protein structure
    view.addModel(exemp_pocket, "pdb") #adds pocket structure
    view.setStyle({"model": 0}, {"cartoon": {"color": "#AB7FE1"}})
    view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "#53C982"}})
    view.zoomTo() #zooms to protein
    showmol(view, height = 200, width = 200)

#example visualization of ligand
with col4:
    st.space("medium")
    st.image("app/images/example_compund.png")

#parameter selection
with col6:
    st.space("small")
    st.select_slider('ADME',
                     options=list(range(0,10)),
                     value = 4)
    st.select_slider('Lipinsky Rule of Five',
                     options = list(range(0, 10)),
                     value = 7)



col1, col2, col3, col4, col5, col6, col7 = st.columns([0.5, 3, 1.5, 3, 1.5, 3, 1])
with col2:
    st.write("<p style = 'color: #AB7FE1; text-align: center;'>Pocket locations and residues</p>", 
             unsafe_allow_html = True)

with col4:
    st.write("<p style = 'color: #AB7FE1; text-align: center;'>list of most promising Ligands</p>", 
             unsafe_allow_html = True)
    
with col6:
    st.write("<p style = 'color: #AB7FE1; text-align: center;" \
             "'>Choose the best criteria for the final selection</p>",
             unsafe_allow_html = True)

st.space("medium")
col1, col2, col3 = st.columns([1.5, 1, 1.5])
with col2:
    if st.button("Start the Docking Pipeline"):
        st.switch_page("pages/0_Input_page.py")

st.space("large")

st.markdown(
    """
    <p style="color: grey; text-align: center;">
        OneDock was created by Thaddeus Kühn, Tine Limberg, Manon Mandernach 
        and Sylviane Verschaeve as part of the Meet-EU project of 2025/26.
        For more information
        <a href="https://github.com/meet-eu-25-26/Heidelberg_Team_2" target="_blank">check out our GitHub page</a>.
    </p>
    """,
    unsafe_allow_html=True
)