#!/usr/bin/env bash
set -e

# --- STEP 1: Install Miniconda (ARM64) if not present ---
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
fi

# Ensure Conda is active in this shell
export PATH="$HOME/miniconda/bin:$PATH"
eval "$($HOME/miniconda/bin/conda shell.bash hook)"

# --- STEP 2: Create or update the minimal conda environment ---
conda env update -f environment-minimal.yml --prune

# --- STEP 3: Activate the environment ---
conda activate version-check-minimal

# --- STEP 4: Run the minimal version-check pipeline ---
nextflow run ../../shared/nextflow/workflows/version-check.nf --outdir results 