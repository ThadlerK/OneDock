FROM continuumio/miniconda3:latest

# Get common tools needed to download files AND install tmux
RUN apt-get update && apt-get install -y \
    wget \
    tmux && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Download the pre-compiled executable, rename it, and make it runnable
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/bin/vina && \
    chmod +x /usr/bin/vina

# Create a dedicated conda environment for your Python tools.
RUN conda create -n adtools -c conda-forge -c bioconda -y \
    python=3.10 \
    meeko \
    rdkit \
    numpy \
    scipy \
    gemmi \
    prody \
    openbabel \
    pandas \
    seaborn \
    matplotlib \
    fpocket \
    pip && \
    conda run -n adtools pip install biobb_structure_utils biobb_vs


# 3. Clean up the conda cache
RUN conda clean -a -y

# Now, 'vina' is "installed" and can be run from anywhere
CMD ["vina", "--help"]