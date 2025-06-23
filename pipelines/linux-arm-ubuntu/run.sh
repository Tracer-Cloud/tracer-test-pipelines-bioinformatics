#!/bin/bash
set -e

echo "[INFO] Starting environment setup..."

# Step 1: Install Java (OpenJDK 17)
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java is already installed."
fi

# Step 2: Install Miniconda (ARM-compatible)
CONDA_DIR="$HOME/miniconda"
if [ ! -d "$CONDA_DIR" ]; then
    echo "[INFO] Installing Miniconda (for ARM)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$CONDA_DIR"
else
    echo "[INFO] Miniconda already installed at $CONDA_DIR"
fi

# Ensure conda is in PATH
export PATH="$CONDA_DIR/bin:$PATH"
if ! command -v conda &> /dev/null; then
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
fi

# Step 3: Install Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mkdir -p "$HOME/.local/bin"
    mv nextflow "$HOME/.local/bin/"
    export PATH="$HOME/.local/bin:$PATH"
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "[INFO] Nextflow is already installed."
fi

# Confirm installs
echo "[✅] Setup complete."
java -version
conda --version
nextflow -version

echo "You can now run:"
echo "   nextflow run main.nf -c nextflow.config"
