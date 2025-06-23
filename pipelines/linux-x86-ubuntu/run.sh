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

# Step 2: Install Miniconda for x86_64
if ! command -v conda &>/dev/null; then
    echo "[INFO] Installing Miniconda for x86_64..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
    export PATH="$HOME/miniconda/bin:$PATH"
    source "$HOME/.bashrc"
else
    echo "[INFO] Conda already installed."
fi

# Step 3: Install Nextflow
if ! command -v nextflow &>/dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

echo "[✅] Setup complete."
echo "Java version:"
java -version
echo "Conda version:"
conda --version || echo "⚠️ Conda not found"
echo "Nextflow version:"
nextflow -version || echo "⚠️ Nextflow not found"
