FROM continuumio/miniconda3:latest

# 1. Install System Tools
RUN apt-get update && apt-get install -y \
    git \
    wget \
    tmux \
    build-essential \
    patch \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Install AutoDock Vina and P2Rank
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    -O /usr/bin/vina && chmod +x /usr/bin/vina

RUN wget https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz \
    -O /tmp/p2rank.tar.gz && \
    mkdir -p /opt && \
    tar -xzf /tmp/p2rank.tar.gz -C /tmp && \
    mv /tmp/p2rank_2.5.1 /opt/ && \
    chmod +x /opt/p2rank_2.5.1/prank && \
    ln -sf /opt/p2rank_2.5.1/prank /usr/bin/p2rank && \
    rm -rf /tmp/p2rank.tar.gz

# 3. Install Mamba explicitly
RUN conda install -n base -c conda-forge mamba -y

# 4. Create Environment & Install Packages in stages
# STAGE A: Create env and install the HEAVY dependencies (Physics/Chem/CUDA)
RUN mamba create -n adtools -c conda-forge -c bioconda --override-channels -y \
    python=3.10 \
    cuda-version=11.8 \
    ambertools=23 \
    openmm \
    openmmforcefields \
    pdbfixer \
    meeko \
    openbabel \
    gemmi \
    biobb_common=5.1.0 \
    biobb_structure_utils=5.1.0

# STAGE B: Install Data Science & Bioinformatics tools into the SAME environment
RUN mamba install -n adtools -c conda-forge -c bioconda --override-channels -y \
    rdkit \
    mdanalysis \
    openjdk=17 \
    numpy \
    scipy \
    prody \
    pandas \
    seaborn \
    matplotlib \
    fpocket \
    snakemake \
    streamlit \
    streamlit-aggrid \
    biopython \
    pip \
    && mamba clean -afy

# 5. Install pip packages (Activating environment via run shell)
SHELL ["/bin/bash", "-c"]
RUN source activate adtools && pip install \
    biobb_vs \
    py3Dmol \
    ipython_genutils \
    stmol \
    bioemu \
    torch \
    posebusters

# 6. Setup Environment Variables
ENV PATH="/opt/conda/envs/adtools/bin:$PATH"
ENV AMBERHOME="/opt/conda/envs/adtools"

# Fix for VS Code Git Issue
RUN git config --global --add safe.directory '*'

# Fix for VS Code Terminal
RUN echo "source activate adtools" >> ~/.bashrc 

# 7. Setup app
WORKDIR /app
COPY . .

# 8. Final config
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility
ENV LD_LIBRARY_PATH="/usr/local/nvidia/lib:/usr/local/nvidia/lib64:$LD_LIBRARY_PATH"
EXPOSE 8501
CMD ["streamlit", "run", "app/Home.py", "--server.port=8501", "--server.address=0.0.0.0"]