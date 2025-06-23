#!/bin/bash
set -e

echo "[INFO] Starting minimal setup for ARM-based Nextflow pipeline..."

# Step 1: Install Java (OpenJDK 17)
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
fi

# Step 2: Install Miniforge (ARM-friendly Conda)
if ! command -v conda &> /dev/null; then
    echo "[INFO] Installing Miniforge (ARM-compatible Conda)..."
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O miniforge.sh
    bash miniforge.sh -b -p "$HOME/miniforge"
    export PATH="$HOME/miniforge/bin:$PATH"
    echo 'export PATH="$HOME/miniforge/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "[INFO] Conda is already installed."
fi

# Step 3: Install Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow is already installed."
fi

# Done
echo "[✅] Environment setup complete."
echo "You can now run:"
echo "   nextflow run main.nf -c nextflow.config"
