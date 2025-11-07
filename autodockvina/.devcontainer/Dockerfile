# Start from a clean miniconda base
FROM continuumio/miniconda3:latest

# Set non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

# Update conda, set up channels, and install autodock-vina
# This single 'autodock-vina' package from conda-forge includes
# the executables (vina, vina_split) AND the python library.
RUN conda update -n base -c defaults conda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda install -y autodock-vina && \
    conda clean -a -y

# Optional: Set a default command to show Vina's help
CMD ["vina", "--help"]