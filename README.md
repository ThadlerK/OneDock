# OneDock Virtual Screening Pipeline

A comprehensive Snakemake + Streamlit pipeline for virtual screening with integrated ADME analysis.

## Features

- **Structure Preparation**: Upload PDB structures or generate with BioEmu
- **Pocket Prediction**: fpocket and P2Rank integration
- **Molecular Docking**: AutoDock Vina with residue-based pocket definition
- **ADME Screening**: SwissADME integration for drug-likeness evaluation
- **Interactive UI**: Streamlit-based web interface

## Project Structure

```
OneDock-Pipeline/
├── config/
│   ├── config.yaml          # Configuration for Snakemake
│   └── samples.csv          # Input metadata (optional)
├── data/
│   ├── inputs/              # Input files (PDB, SMILES)
│   │   └── library_split/   # Individual ligand files
│   ├── interim/             # Intermediate files (PDBQT conversions)
│   └── results/             # Output files (docking scores, ADME)
├── workflow/
│   ├── Snakefile            # Main Snakemake workflow
│   ├── envs/                # Conda environments
│   ├── rules/               # Modularized rules
│   └── scripts/             # Pipeline scripts
│       ├── bioemu_pipeline.py
│       ├── convert_smi_to_pdbqt.py
│       └── docking_script_with_residues.sh
├── app/
│   ├── Home.py              # Main Streamlit app entry point
│   ├── utils.py             # Helper functions
│   ├── pages/               # Multi-page layout
│   │   ├── 1_Structure_Preparation.py
│   │   ├── 2_Pocket_Prediction.py
│   │   ├── 3_Docking.py
│   │   ├── 4_Results.py
│   │   └── 5_ADME_Screening.py
│   └── tools/               # Analysis tools
│       ├── fpocket.py
│       ├── P2Rank_filtering.py
│       ├── P2Rank_to_PDB.py
│       ├── pocket_comparison.py
│       └── swissadme_client.py
├── docs/
│   └── ADME_INTEGRATION.md  # ADME screening documentation
├── .gitignore
├── README.md
└── requirements.txt         # Python dependencies
```

## Quick Start

### 1. Start the Application
```bash
streamlit run app/Home.py --server.address=0.0.0.0 --server.port=8501
```

### 2. Pipeline Workflow

1. **Input & Visualization** (Home.py)
   - Upload receptor PDB or generate with BioEmu
   - Upload ligand library (SMILES format)

2. **Structure Preparation** (Page 1)
   - Convert receptor to PDBQT format
   - Prepare ligands for docking

3. **Pocket Prediction** (Page 2)
   - Run fpocket or P2Rank
   - Visualize predicted binding sites

4. **Docking** (Page 3)
   - Configure docking parameters
   - Launch AutoDock Vina pipeline
   - Define pocket by residues

5. **Results** (Page 4)
   - View docking scores
   - Rank ligands by binding affinity

6. **ADME Screening** (Page 5) 🆕
   - Select top candidates
   - Submit to SwissADME
   - Filter by drug-likeness rules
   - Combined docking + ADME analysis

## ADME Integration

The pipeline now includes comprehensive ADME property evaluation:

- **Drug-likeness**: Lipinski's Rule of Five, Veber rules, Ghose filter
- **Pharmacokinetics**: GI absorption, BBB permeability, bioavailability
- **Filtering**: Remove compounds with PAINS alerts
- **Visualization**: Interactive plots for property analysis

See [ADME_INTEGRATION.md](docs/ADME_INTEGRATION.md) for detailed documentation.

## Requirements

- Python 3.10+
- Streamlit
- Snakemake
- AutoDock Vina
- Open Babel
- fpocket / P2Rank
- Pandas, Plotly

## Installation

```bash
pip install -r requirements.txt
```

## Configuration

Edit `config/config.yaml` to customize:
- Receptor path
- Pocket residues
- Grid size and exhaustiveness
- Docking engine selection

## Output Files

```
data/results/
├── docking_report.csv              # Docking scores
├── swissadme_results.csv           # ADME properties
├── combined_docking_adme.csv       # Integrated analysis
├── poses/                          # Docked conformations
├── logs/                           # Vina logs
└── stats/                          # Per-ligand statistics
```

## Citation

If you use this pipeline, please cite:
- SwissADME: Daina et al. (2017) Scientific Reports 7:42717
- AutoDock Vina: Trott & Olson (2010) J Comput Chem 31:455-461

## License

MIT License