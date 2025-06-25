#!/usr/bin/env bash
set -e

if [ ! -d "$HOME/miniconda" ]; then
    echo "[INFO] Installing Miniconda (ARM64)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"
    if ! grep -q 'miniconda/bin' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
    fi
    "$HOME/miniconda/bin/conda" init bash
    echo "[INFO] Miniconda installed. Please run 'source ~/.bashrc' or restart your shell to use conda."
fi

export PATH="$HOME/miniconda/bin:$PATH"
eval "$("$HOME/miniconda/bin/conda" shell.bash hook)"

if conda env list | grep -q "^linux-arm-ubuntu-minimal[[:space:]]"; then
    echo "[INFO] Conda environment 'linux-x86-ubuntu-minimal' already exists. Using it as is."
else
    echo "[INFO] Creating conda environment 'linux-x86-ubuntu-minimal'..."
    conda env update -f environment-minimal.yml
fi

conda activate linux-x86-ubuntu-minimal

nextflow -log logs/nextflow.log run main.nf --outdir results
