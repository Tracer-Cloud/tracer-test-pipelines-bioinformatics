#!/bin/bash

set -e

echo "[INFO] Starting pipeline setup..."

# Step 1: Install Conda (if not available)
if ! command -v conda &> /dev/null; then
    echo "[INFO] Conda not found. Installing Miniconda for ARM..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "[INFO] Conda is already installed."
fi

# Step 2: Initialize Conda (just in case)
eval "$($HOME/miniconda/bin/conda shell.bash hook)" || true

# Step 3: Run the pipeline
echo "[INFO] Running Nextflow pipeline..."
nextflow run main.nf -c nextflow.config
