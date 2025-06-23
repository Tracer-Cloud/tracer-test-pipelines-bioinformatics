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

# Step 2: Install Miniforge (ARM-compatible Conda)
if ! command -v conda &> /dev/null; then
    echo "[INFO] Installing Miniforge (ARM-compatible Conda)..."
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O miniforge.sh
    bash miniforge.sh -b -u -p "$HOME/miniforge"
    export PATH="$HOME/miniforge/bin:$PATH"
    echo 'export PATH="$HOME/miniforge/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "[INFO] Miniforge already exists at $HOME/miniforge"
    export PATH="$HOME/miniforge/bin:$PATH"
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

# Final Check
echo "[✅] Environment setup complete."
echo "Java version:"
java -version
echo "Nextflow version:"
nextflow -version
echo "You can now run:"
echo "   nextflow run main.nf -c nextflow.config"
