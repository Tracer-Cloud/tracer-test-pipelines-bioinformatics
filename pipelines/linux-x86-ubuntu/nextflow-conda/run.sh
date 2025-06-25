#!/bin/bash
set -e

echo "[INFO] Minimal setup for Ubuntu x86_64..."

# --- STEP 1: Install Java ---
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java already installed."
fi

# --- STEP 2: Install Miniconda (x86_64) ---
if [ ! -d "$HOME/miniconda" ]; then
    echo "[INFO] Installing Miniconda (x86_64)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"

    # Add conda to PATH in .bashrc if not already present
    if ! grep -q 'miniconda/bin' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
    fi

    # Initialize conda for bash shell
    "$HOME/miniconda/bin/conda" init bash

    echo "[INFO] Miniconda installed. Please run 'source ~/.bashrc' or restart your shell to use conda."
fi

# Ensure Conda is active in this shell
export PATH="$HOME/miniconda/bin:$PATH"
eval "$($HOME/miniconda/bin/conda shell.bash hook)"

# --- STEP 3: Create conda environment with bioinformatics tools ---
if ! conda env list | grep -q "rnaseq-tools"; then
    echo "[INFO] Creating conda environment with bioinformatics tools..."
    conda env create -f environment.yml
    echo "[INFO] Conda environment 'rnaseq-tools' created successfully."
    echo "[INFO] To activate: conda activate rnaseq-tools"
else
    echo "[INFO] Conda environment 'rnaseq-tools' already exists."
fi

# --- STEP 4: Install Nextflow ---
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

# --- Final check ---
echo "[✅] Environment ready:"
echo
java -version
conda --version
nextflow -version

echo ""
echo "[INFO] Setup complete! To run nf-core pipelines:"
echo "1. Activate the conda environment: conda activate rnaseq-tools"
echo "2. Run your pipeline: nextflow run nf-core/rnaseq -c custom.config -profile test --outdir results"
echo ""
echo "[INFO] To use conda in new terminal sessions, run:"
echo "source ~/.bashrc"
echo ""
echo "Or restart your terminal/shell session."
