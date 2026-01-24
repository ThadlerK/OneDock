FROM continuumio/miniconda3:latest

# 1. Install System Tools
RUN apt-get update && apt-get install -y \
    git \
    wget \
    tmux \
    build-essential \
    openjdk-17-jre-headless \
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

# 3. Install Libmamba Solver explicitly
RUN conda install -n base -c conda-forge conda-libmamba-solver -y && \
    conda config --set solver libmamba

# 4. Create the Conda environment (FORCE libmamba here)
# We group the pip installs here to keep the layer size smaller
RUN conda create -n adtools -c conda-forge -c bioconda --solver=libmamba -y \
    python=3.10 \
    cuda-version=11.8 \
    ambertools=23 \
    openmmforcefields \
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
    biobb_common=5.1.0 \
    biobb_structure_utils=5.1.0 \
    snakemake \
    streamlit \
    streamlit-aggrid \
    biopython \
    openmm \
    pdbfixer \
    pip \
    && conda clean -afy

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
EXPOSE 8501
CMD ["streamlit", "run", "app/Home.py", "--server.port=8501", "--server.address=0.0.0.0"]