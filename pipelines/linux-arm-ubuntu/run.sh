#!/bin/bash
set -e

echo "[INFO] Starting minimal setup for ARM-based Nextflow pipeline..."

# Step 1: Install Java (OpenJDK 17)
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java is already installed."
fi

# Step 2: Install Miniforge (ARM-compatible Conda)
if ! command -v conda &> /dev/null; then
    if [ ! -d "$HOME/miniforge" ]; then
        echo "[INFO] Installing Miniforge (ARM-compatible Conda)..."
        wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O miniforge.sh
        bash miniforge.sh -b -p "$HOME/miniforge"
    else
        echo "[INFO] Miniforge directory already exists, skipping install."
    fi

    export PATH="$HOME/miniforge/bin:$PATH"
    echo 'export PATH="$HOME/miniforge/bin:$PATH"' >> ~/.bashrc
else
    echo "[INFO] Conda (Miniforge) is already available."
fi

# Step 3: Install Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow "$HOME/.local/bin/"

    # Ensure ~/.local/bin is in PATH
    mkdir -p "$HOME/.local/bin"
    export PATH="$HOME/.local/bin:$PATH"
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
else
    echo "[INFO] Nextflow is already installed."
fi

# Final check
echo
echo "[✅] Environment setup complete."
echo "Java version: $(java -version 2>&1 | head -n 1)"
echo "Nextflow version: $(nextflow -version | head -n 1)"
echo
echo "You can now run the pipeline with:"
echo "   nextflow run main.nf -c nextflow.config"
