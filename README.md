![Alt text](OneDock/app/images/logo_dunkel.png)
<p align = "center"><strong>The easiest way to find the most promising Ligands for your Protein</strong></p>

Welcome to your protein-ligand matchmaker!
If you want to find out which ligands might bind to your protein, you have come to the right place. OneDock provides you with a variety of possibilities to screen your library of small molecules for the best candidates for your protein of interest, even if you don't know the binding pockets of your protein or its experimental structure.



<table align = "center" style = "border-collapse: collapse;">
  <tr>
    <td align = "center">
      <img src = "OneDock/app/images/specificity.svg" ><br>
      <sub>Ligands with the highest <strong>specificity</strong> and <strong>selectivity</strong></sub>
    </td>
    <td align = "center">
      <img src = "OneDock/app/images/unknown_sites.svg"><br>
      <sub>Identify <strong>unknown binding sites</strong></sub>
    </td>
    <td align = "center">
      <br>
      <br>
      <br>
      <img src = "OneDock/app/images/protein_sequence.svg"><br>
      <br>
      <br>
      <sub>Work with protein <strong>sequence data</strong></sub>
    </td>
  </tr>
</table>



This Readme will take you through the workings of the pipeline. We'll show you what possibilities you have, how to use it on your own device, and what you need to consider when using the pipeline.
<br>
<br>

## How to find your way through the Pipeline 
![Alt text](OneDock/app/images/pipeline.svg)

As you can see, there are a lot of different possibilities, so let's walk through them. 

 
<details>
  <summary> Scenario 1:You already know the protein structure and binding pocket</summary>
    If you allready know the protein binding region of your ligand candidate (e.g. a highly functional part of the protein) and only want ligands that bind to this pocket, you can just upload your protein and your ligand library on the first input page and directly proceed to docking. 
</details><br>

<details>
  <summary> Scenario 2: You don't know the binding pocket, but you know the protein structure</summary>
    Then you can try out our pocket detection page. There, you can predict the pockets using detection tools <em>fpocket</em> and <em>P2Rank</em>. Using filters (that you can set yourself) and by looking for overlaps between the two methods, we'll compute the best pockets for you, so that you can proceed with docking. 
</details><br>

<details>
  <summary> Scenario 3: I only have the sequence of my protein </summary>
    We'll use <em>BioEmu</em> to predict the structure and use the pocket detection page. There, you can predict the pockets using detection tools <em>fpocket</em> and <em>P2Rank</em>. Using filters (that you can set) and by looking for overlaps between the two methods, we'll compute the best pockets for you, so that you can proceed with docking. 
</details><br>


On the docking page, you will use <em>AutoDock Vina</em> to dock the ligands to the (pre)defined pocket. It will give out the best ligand candidates which you can filter on the ranking site using <em>ADME</em> and <em>Lipinski's rule of five</em>. 

In the end you will be provided with your pocket residues, your final ligand candidates and their respective filtering scores. Additionally, OneDock will visualize the ligands candidates and how they bind the protein. 
By looking at the different scores you can then identify the ligand that fits your requirements best. 

<details>
  <summary> What about other targets?</summary>
    You can upload an additional reference strucuture on the upload page. OneDock will repeat the procedure for this protein and compare the docking results to the results of your protein of interest. This way you get the ligands that are specific for your protein. 
</details><bp><br>

## How is the code structured and how can I use it on my device?
### General structure
This is the structure of OneDock and this GitHub page:<bp>
```text
OneDock
├── .devcontainer/          
│   ├── devcontainer.json
│   └── Dockerfile
├── app/
│   ├── .streamlit/         # streamlit layout
│   │   └── config.toml
│   ├── images/             #images needed for the app
│   ├── pages/              # streamlit pages
│   │   ├── 1_Structure_Preparation.py
│   │   ├── 2_Pocket_Prediction.py
│   │   ├── 3_Docking.py
│   │   ├── 4_Results.py
│   │   └── 5_ADME_Screening.py
│   ├── Home.py              # Main streamlit app
│   └── utils.py             # Helper functions 
├── config/ 
│   └── config.yaml          # Configuration for Snakemake
├── data/
│   ├── inputs/              # Input files
│   ├── interim/             # Interim files 
│   └── results/             # Output files 
├── workflow/
│   ├── scripts/              
│   │   ├── bioemu_pipeline.py
│   │   ├── convert_smi_to_bdbqt.py
│   │   ├── docking_script_with_residues.py
│   │   ├── fix_atom_names.py
│   │   ├── mmpbsa_worker.sh
│   │   ├── run_openmm_md_old.py
│   │   └── run_openmm_md.py
│   └── Snakefile            # Main entry point for Snakemake
└── .gitignore
README.md
images_readme/              

```
<br>

### How to use OneDock
For the ideal use of OneDock install [VS Code](https://code.visualstudio.com/download) and [docker](https://www.docker.com/). You can then clone the repository:

```bash
git clone https://github.com/meet-eu-25-26/Heidelberg_Team_2.git
```
Open the OneDock folder in VS code and choose *Dev Containers: Rebuild and Reopen in Container*.
Once you are in your container you can start OneDock:
```bash
streamlit run app/Home.py
```
OneDock will guide you through the steps according to your wishes. 
You will need GPU access to use the BioEmu and MM/GBSA tools.<br>



### The Tools

#### **Streamlit**
<p align = "center">
<img src = "images_readme/streamlit-logo.png" width="20%"></p>

We are using [streamlit](https://streamlit.io/) to build OneDock. It is a tool that helps you create an interactive app that allows you to set parameters for computation and visualize your results in different ways. This way OneDock is as user-friendly as possible, to provide you with the best experience possible!

#### **Snakemake** 
<p align = "center">
<img src = "images_readme/snakemake.jpg" width = "30%"></p>

The pipeline is organized in [snakemake](https://snakemake.github.io/). It allows us to create a workflow system in which files are processed by functions (so called rules) while considering their dependencies. This way, we get a safe and reproducible workflow that is easily integrateable from python. 

#### **Docker**
<p align = "center">
<img src = "images_readme/docker.png" width="15%"></p> 

We use a [docker](https://www.docker.com/) container to bundle all the dependencies of OneDock. All of the different libraries, system tools,... necessary to run it are defined in a dockerfile from which an image is built. This way, when you use the pipleine, you won't have to worry about installing everything. 

#### **BioEmu**
<p align = "center">
<img src = "OneDock/app/images/bioemu.svg" width="10%"></p> 
BioEmu is a generative structural modeling approach used to  sample  biomolecular conformations. In this project, BioEmu was used only to generate the  initial protein structures. No dynamics or free energy calculations were performed with BioEmu; it served exclusively as a tool for initial structure sampling.

#### **fpocket**
<p align = "center">
<img src = "OneDock/app/images/fpocket.svg" width="15%"></p> 
fpocket is a pocket detection tool  based on geometric methods such as Voronoi tessellation and α-spheres to identify surface cavities efficiently. It clusters these geometric features into pockets and scores them to estimate properties like druggability and pocket volume <sup>2</sup>.

#### **P2Rank**
<p align = "center">
<img src = "OneDock/app/images/p2rank.svg" width="20%"></p> 
P2Rank is a pocket detection tool that samples points on the protein’s solvent‑accessible surface and assigns each a ligandability score using a random forest model trained on known protein–ligand complexes. Points with high ligandability are clustered into pockets, which are then  scored and ranked by combining the scores of their points <sup>3</sup>.

#### **AutoDock Vina**
<p align = "center">
<img src = "OneDock/app/images/vina.svg" width="7.5%"></p> 
AutoDock Vina is a fast molecular docking tool that predicts binding poses and binding affinities of small molecules to protein targets. It explores a user-defined binding site, generates multiple ligand poses, and scores them using an empirical scoring function <sup>4</sup>.

#### **SwissADME**
<p align = "center">
<img src = "OneDock/app/images/adme.svg" width="20%"></p> 
SwissADME is a web-based tool that uses a variety of predictive models to compute physicochemical descriptors of small molecules from their structure. Thereby it predicts ADME properties, drug-likeness, pharmacokinetics and medicinal chemistry features <sup>5</sup>.

#### **PoseBusters**
<p align = "center">
<img src = "OneDock/app/images/posebusters.svg" width="7.5%"></p> 
PoseBusters is a Python toolkit that validates the physical and chemical plausibility of protein-ligand complexes <sup>7</sup>.

#### **py3Dmol** 
py3Dmol is a python toolkit that enables an interactive 3D visualization of protein-ligand complexes <sup>8</sup>.


#### **OpenMM**
<p align = "center">
<img src = "OneDock/app/images/openmm.svg" width="12%"></p> 
OpenMM is an open-source, high-performance molecular simulation toolkit designed for running molecular dynamics (MD) simulations, with native support for GPU acceleration. It provides implementations of modern force fields and integrators and allows flexible scripting of simulation protocols in Python.
**Setup:** Protein - Ligand complexes were parameterized with the AMBER ff14SB force field for the protein and GAFF2 for the ligand (ligand parameters generated with Antechamber and missing terms via Parmchk2). Molecular dynamics simulations were carried out in OpenMM using an implicit solvent Generalized Born OBC2 model, including energy minimization followed by heating, equilibration, and production runs with gradually released backbone restraints. Trajectories from production were post-processed with MMPBSA.py using the GBSA model to estimate binding energies from snapshots<sup>9</sup>.

#### **MM/GBSA** 
MM/GBSA (Molecular Mechanics / Generalized Born Surface Area) is an end-point free energy method used to estimate binding free energies from molecular dynamics trajectories. It combines molecular mechanics energy terms (electrostatics and van der Waals interactions) with an implicit solvent model based on Generalized Born theory and a nonpolar surface area contribution. MM/GBSA is closely related to MM/PBSA, which instead uses a numerical Poisson–Boltzmann solver for the polar solvation energy. While MM/PBSA can be more detailed, it is also significantly more computationally demanding and sensitive to numerical settings. We chose MM/GBSA because it is faster, more robust, and well suited for screening and ranking many protein-ligand complexes.


<br>
<br>


### References
1. Lewis, S., Hempel, T., Jiménez-Luna, J., Gastegger, M., Xie, Y., Foong, A. Y., ... & Noé, F. (2025). Scalable emulation of protein equilibrium ensembles with generative deep learning. Science, 389(6761), eadv9817.

2. Le Guilloux, V., Schmidtke, P., & Tuffery, P. (2009). Fpocket: an open source platform for ligand pocket detection. BMC bioinformatics, 10(1), 168.

3. Krivák, R., & Hoksza, D. (2018). P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. Journal of cheminformatics, 10(1), 39. 

4. Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. Journal of computational chemistry, 31(2), 455-46

5. Daina, A., Michielin, O., & Zoete, V. (2017). SwissADME: a free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. Scientific reports, 7(1), 42717.

6. Buttenschoen, M., Morris, G.M. & Deane, C.M. 
PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. 
Chem. Sci. 15, 3130-3139 (2024)

7. Rego, N. and Koes, D. 3Dmol.js: molecular visualization with WebGL. 
Bioinformatics 31, 1322-1324 (2015)

8. An Efficient Program for End-State Free Energy Calculations
Bill R. Miller III, T. Dwight McGee Jr., Jason M. Swails, Nadine Homeyer, Holger Gohlke, and Adrian E. Roitberg Journal of Chemical Theory and Computation 2012 8 (9), 3314-3321

9.  P. Eastman, R. Galvelis, R. P. Peláez, C. R. A. Abreu, S. E. Farr, E. Gallicchio, A. Gorenko, M. M. Henry, F. Hu, J. Huang, A. Krämer, J. Michel, J. A. Mitchell, V. S. Pande, J. PGLM Rodrigues, J. Rodriguez-Guerra, A. C. Simmonett, S. Singh, J. Swails, P. Turner, Y. Wang, I. Zhang, J. D. Chodera, G. De Fabritiis, and T. E. Markland. “OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials.” J. Phys. Chem. B 128(1), pp. 109-116 (2023).


<br>
<br>
<br>

OneDock was created by Thaddeus Kühn, Tine Limberg, Manon Mandernach and Sylviane Verschaeve as part of the Meet-EU project of 2025/26. <br>
The project is protected under a MIT license.
