# Use an official Ubuntu base image
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

# Install basic utilities
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    sudo \
    git \
    && apt-get clean

# Install Miniconda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh -O anaconda.sh && \
    bash anaconda.sh -b -p /opt/conda && \
    rm anaconda.sh && \
    /opt/conda/bin/conda init bash

# Install Salmon in a new Conda environment
RUN /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda create -n salmon -y salmon

# Set the entrypoint to bash
ENTRYPOINT ["bash"]
