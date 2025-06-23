#!/bin/bash
set -e

echo "[INFO] Starting environment setup..."

# Detect shell profile (bash assumed)
PROFILE="$HOME/.bashrc"

# Step 1: Install Java (OpenJDK 17)
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java is already installed."
fi

# Step 2: Install Miniconda (ARM)
CONDA_DIR="$HOME/miniconda"
if [ ! -d "$CONDA_DIR" ]; then
    echo "[INFO] Installing Miniconda (for ARM)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$CONDA_DIR"
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$PROFILE"
    export PATH="$HOME/miniconda/bin:$PATH"
else
    echo "[INFO] Miniconda already installed at $CONDA_DIR"
    export PATH="$HOME/miniconda/bin:$PATH"
fi

# Check Conda
if ! command -v conda &> /dev/null; then
    echo "⚠️ Conda not found even after install. Try running: source ~/.bashrc"
else
    echo "[INFO] Conda version: $(conda --version)"
fi

# Step 3: Install Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mkdir -p "$HOME/.local/bin"
    mv nextflow "$HOME/.local/bin/"
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$PROFILE"
    export PATH="$HOME/.local/bin:$PATH"
else
    echo "[INFO] Nextflow already installed: $(nextflow -version | head -n 1)"
fi

# Final echo
echo -e "\n[✅] Setup complete.\n"
java -version
conda --version || echo "⚠️ Conda not available in PATH"
nextflow -version || echo "⚠️ Nextflow not available in PATH"

echo -e "\nYou can now run:\n  nextflow run main.nf -c nextflow.config"
