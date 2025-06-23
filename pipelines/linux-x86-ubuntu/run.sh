#!/bin/bash
set -e

echo "[INFO] Starting environment setup for linux-x86-ubuntu..."

# Step 1: Install Java
if ! command -v java &>/dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update
    sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java is already installed."
fi

# Step 2: Install Miniconda (clean install)
CONDA_DIR="$HOME/miniconda"

if [ ! -x "$CONDA_DIR/bin/conda" ]; then
    echo "[INFO] Installing Miniconda for x86_64..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    rm -rf "$CONDA_DIR"
    bash miniconda.sh -b -p "$CONDA_DIR"
    echo "export PATH=\"$CONDA_DIR/bin:\$PATH\"" >> ~/.bashrc
    export PATH="$CONDA_DIR/bin:$PATH"
    source ~/.bashrc
else
    echo "[INFO] Conda already installed."
    export PATH="$CONDA_DIR/bin:$PATH"
fi

# Step 3: Install Nextflow
if ! command -v nextflow &>/dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

echo "[✅] Setup complete."
echo
echo "🔍 Versions:"
java -version
conda --version || echo "⚠️ Conda not found in PATH"
nextflow -version || echo "⚠️ Nextflow not found"
