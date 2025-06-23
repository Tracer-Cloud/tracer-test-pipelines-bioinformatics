#!/bin/bash
set -e

echo "[INFO] Starting environment setup..."

# Step 1: Install Java
if ! command -v java &>/dev/null; then
    echo "[INFO] Installing OpenJDK..."
    sudo apt-get update
    sudo apt-get install -y default-jdk
else
    echo "[INFO] Java is already installed."
fi

# Step 2: Install Miniconda (x86_64)
if ! command -v conda &>/dev/null; then
    echo "[INFO] Installing Miniconda for x86_64..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"
    export PATH="$HOME/miniconda/bin:$PATH"
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "[INFO] Conda already installed."
fi

# Step 3: Install Nextflow (clean and force it!)
if ! command -v nextflow &>/dev/null; then
    echo "[INFO] Installing Nextflow..."
