#!/bin/bash
set -e

echo "[INFO] Setting up Nextflow pipeline with Conda for Ubuntu ARM..."

# --- STEP 1: Install Docker ---
if ! command -v docker &> /dev/null; then
    echo "[INFO] Installing Docker..."
    sudo apt-get update
    sudo apt-get install -y docker.io
    sudo systemctl start docker
    sudo systemctl enable docker
    sudo usermod -aG docker $USER
    echo "[INFO] Docker installed. You may need to log out and back in for group changes to take effect."
else
    echo "[INFO] Docker already installed."
fi

# --- STEP 2: Install Java ---
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java already installed."
fi

# --- STEP 3: Install Miniconda (ARM64) ---
if [ ! -d "$HOME/miniconda" ]; then
    echo "[INFO] Installing Miniconda (ARM64)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"

    # Add conda to PATH in .bashrc if not already present
    if ! grep -q 'miniconda/bin' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
    fi

    "$HOME/miniconda/bin/conda" init bash

    echo "[INFO] Miniconda installed. Please run 'source ~/.bashrc' or restart your shell to use conda."
else
    echo "[INFO] Miniconda already installed."
fi

# Ensure Conda is active in this shell
export PATH="$HOME/miniconda/bin:$PATH"
eval "$($HOME/miniconda/bin/conda shell.bash hook)"

# --- STEP 4: Create conda environment with all dependencies ---
echo "[INFO] Creating conda environment 'nextflow-pipeline' with all dependencies..."
conda env create -f environment-minimal.yml

# --- STEP 5: Install Nextflow in the conda environment ---
echo "[INFO] Installing Nextflow in the conda environment..."
conda activate rnaseq-minimal
conda install -c bioconda nextflow -y

echo "[INFO] Environment ready!"
echo
echo "[INFO] ========================================="
echo "[INFO] SETUP COMPLETE!"
echo "[INFO] ========================================="
echo
echo "[INFO] To use the pipeline:"
echo
echo "1. Activate the conda environment:"
echo "   conda activate nextflow-pipeline"
echo
echo "2. Verify the environment:"
echo "   conda list"
echo
echo "3. Run the pipeline:"
echo "   nextflow run main.nf --outdir results"
echo
echo "4. Run nf-core RNA-seq pipeline:"
echo "   nextflow run nf-core/rnaseq -r 3.14.0 -c custom.config -profile test --outdir results -resume"
echo
echo "[INFO] To use conda in new terminal sessions, run:"
echo "source ~/.bashrc"
echo
echo "Or restart your terminal/shell session."
echo
echo "[INFO] =========================================" 