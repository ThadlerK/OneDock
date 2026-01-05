FROM continuumio/miniconda3:latest

# 1. Install System Tools
RUN apt-get update && apt-get install -y \
    git \
    wget \
    tmux \
    build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 2. Install AutoDock Vina (Binary) and P2Rank
# We put it in /usr/bin so it's globally available
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/bin/vina && \
    chmod +x /usr/bin/vina

RUN wget https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz \
        -O /tmp/p2rank.tar.gz && \
    tar -xzf /tmp/p2rank.tar.gz -C /opt && \
    ln -s /opt/p2rank_2.5.1/p2rank.sh /usr/bin/p2rank && \
    rm /tmp/p2rank.tar.gz

# 3. Create the Conda Environment
# We install Streamlit and Snakemake HERE so they can see the bio-tools
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
    biobb_structure_utils \
    snakemake \
    streamlit \
    pip 

# 4. Install Pip Packages into the 'adtools' environment
# Using 'conda run' ensures they go into the right place
RUN conda run -n adtools pip install \
    biobb_vs \
    py3Dmol \
    ipython_genutils \
    stmol

# 5. [CRITICAL] Make 'adtools' the default environment
# By adding the environment's bin folder to PATH, we "activate" it permanently.
ENV PATH="/opt/conda/envs/adtools/bin:$PATH"

# 6. Setup App
WORKDIR /app
COPY . .

# 7. Final Config
EXPOSE 8501

CMD ["bash"]
# Correct CMD: Just run Streamlit
#CMD ["streamlit", "run", "app/Home.py", "--server.address=0.0.0.0"]