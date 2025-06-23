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
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
fi

# Ensure Conda is active in this shell
export PATH="$HOME/miniconda/bin:$PATH"
eval "$($HOME/miniconda/bin/conda shell.bash hook)"

# --- STEP 3: Install Nextflow ---
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
