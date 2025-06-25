#!/bin/bash
# Source this script to activate conda in your current shell session
# Usage: source activate-conda.sh

if [ -d "$HOME/miniconda" ]; then
    export PATH="$HOME/miniconda/bin:$PATH"
    eval "$($HOME/miniconda/bin/conda shell.bash hook)"
    echo "[INFO] Conda activated in current shell session."
    echo "Conda version: $(conda --version)"
    echo ""
    echo "You can now run: nextflow run main.nf --outdir results"
else
    echo "[ERROR] Miniconda not found at $HOME/miniconda"
    echo "Please run ./run.sh first to install Miniconda."
fi
