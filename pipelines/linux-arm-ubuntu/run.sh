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

# Step 2: Install Miniforge (Conda for ARM)
if ! command -v conda &> /dev/null; then
    if [ ! -d "$HOME/miniforge" ]; then
        echo "[INFO] Installing Miniforge..."
        wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O miniforge.sh
        bash miniforge.sh -b -p "$HOME/miniforge"
    else
        echo "[INFO] Miniforge already exists at $HOME/miniforge"
    fi

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
    mv nextflow "$HOME/.local/bin/"

    # Ensure it's in PATH
    mkdir -p "$HOME/.local/bin"
    export PATH="$HOME/.local/bin:$PATH"
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "[INFO] Nextflow is already installed."
fi

# Confirm everything is ready
echo "[✅] Environment setup complete."
echo
java -version
conda --version
nextflow -version

echo
echo "[ℹ️] You can now run:"
echo "     nextflow run main.nf -c nextflow.config"
