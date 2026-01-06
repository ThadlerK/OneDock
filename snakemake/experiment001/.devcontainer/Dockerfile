FROM continuumio/miniconda3:latest

# 1. Install System Tools
RUN apt-get update && apt-get install -y \
    git \
    wget \
    tmux \
    build-essential \
    openjdk-17-jre-headless \
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


# 3. Make conda non-interactive and up-to-date
RUN conda update -n base -c defaults -y conda && \
    conda config --set always_yes yes --set changeps1 no

# 4. Create the Conda environment
RUN conda create -n adtools -c conda-forge -c bioconda \
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
    biobb_structure_utils \
    snakemake \
    streamlit \
    biopython \
    pip

# 5. Install pip packages
RUN /opt/conda/envs/adtools/bin/pip install \
    biobb_vs \
    py3Dmol \
    ipython_genutils \
    stmol

# 6. Clean conda caches (IMPORTANT)
RUN conda clean -afy

# 7. Activate adtools by default
ENV PATH="/opt/conda/envs/adtools/bin:$PATH"

# 8. Setup app
WORKDIR /app
COPY . .

# 9. Final config
EXPOSE 8501
CMD
