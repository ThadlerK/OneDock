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
col1, col2, col3 = st.columns (3)
with col1:
    with st.container(border = True, horizontal_alignment = "center", height = 240):
        st.image("app/images/specificity.svg")
        st.write("<p style = 'color: #AB7FE1; text-align: center;'>Find the ligands with the highest specificity for your protein</p>", 
                 unsafe_allow_html = True)
with col2:
    with st.container(border = True, horizontal_alignment = "center", height = 240):
        st.image("app/images/unknown_sites.svg")
        st.write("<p style = 'color: #AB7FE1; text-align: center;'>Identify unknown binding sites</p>",
                 unsafe_allow_html = True)
with col3:
    with st.container(border = True, horizontal_alignment = "center", height = 240):
        st.space("medium")
        st.image("app/images/protein_sequence.svg")
        st.space("medium")
        st.write("<p style = 'color: #AB7FE1; text-align: center;'>Work with proteins without a defined structure</p>",
                 unsafe_allow_html = True)

st.space('medium')

st.write('To do this you will be guided through the following pipeline:')

st.image('app/images/pipeline.png')

st.space('large')

st.write('Once you have completed the pipeline OneDock will provide you with \n' \
         'extensive details about your ligands:')

col1, col2, col3 = st.columns(3)

#example visualization of pockets
with col1:
    with st.container(border = True, horizontal_alignment = "center", height = 240):
        with open("app/images/example_protein.pdb", "r") as f:
            exemp_protein = f.read()
        with open("app/images/example_pocket.pdb", "r") as f:
            exemp_pocket = f.read()
        view = py3Dmol.view(width=150, height= 150) #sets the room and its size
        view.addModel(exemp_protein, "pdb") #adds protein structure
        view.addModel(exemp_pocket, "pdb") #adds pocket structure
        view.setStyle({"model": 0}, {"cartoon": {"color": "#AB7FE1"}})
        view.setStyle({"model": 1}, {"sphere": {"radius": 1.0, "color": "#53C982"}})
        view.zoomTo() #zooms to protein
        showmol(view, height = 150, width = 150)
        st.write("<p style = 'color: #AB7FE1; text-align: center;'>Pocket locations and residues</p>", 
                 unsafe_allow_html = True)

#example visualization of ligand
with col2:
    with st.container(border = True, horizontal_alignment = "center", height = 240):
        st.space("small")
        st.image("app/images/example_compound.png")
        st.space("small")
        st.write("<p style = 'color: #AB7FE1; text-align: center;'>List of the most promising ligands</p>", 
                 unsafe_allow_html = True)

#parameter selection
with col3:
    with st.container(border = True, horizontal_alignment = "center", height = 240):
        st.space("small")
        st.image("app/images/ADME_regler.png")
        st.write("<p style = 'color: #AB7FE1; text-align: center;" \
                 "'>Choose the best criteria for the final selection</p>",
                 unsafe_allow_html = True)



st.space("medium")
col1, col2, col3 = st.columns([1.5, 1, 1.5])
with col2:
    if st.button("Start OneDock"):
        st.switch_page("pages/0_Input.py")

st.space("large")

st.markdown(
    """
    <p style="color: grey; text-align: center;">
        OneDock was created by Thaddeus Kühn, Tine Limberg, Manon Mandernach 
        and Sylviane Verschaeve as part of the Meet-EU project of 2025/26. <br>
        For more information
        <a href="https://github.com/meet-eu-25-26/Heidelberg_Team_2" target="_blank">check out our GitHub page</a>.
    </p>
    """,
    unsafe_allow_html=True
)