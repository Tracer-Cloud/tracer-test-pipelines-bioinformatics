#!/bin/bash
set -e

echo "[INFO] Minimal setup for Ubuntu x86_64..."

# --- STEP 1: Install Docker ---
if ! command -v docker &> /dev/null; then
    echo "[INFO] Installing Docker..."
    sudo apt-get update
    sudo apt-get install -y docker.io
    sudo systemctl start docker
    sudo systemctl enable docker
    sudo usermod -aG docker $USER
    echo "[INFO] Docker installed. You may need to log out and back in for group changes to take effect."
else
    echo "[INFO] Docker already installed."
fi

# --- STEP 2: Install Java ---
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java already installed."
fi

# --- STEP 3: Install Miniconda (x86_64) ---
if [ ! -d "$HOME/miniconda" ]; then
    echo "[INFO] Installing Miniconda (x86_64)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"

    # Add conda to PATH in .bashrc if not already present
    if ! grep -q 'miniconda/bin' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
    fi

    # Initialize conda for bash shell
    "$HOME/miniconda/bin/conda" init bash

    echo "[INFO] Miniconda installed. Please run 'source ~/.bashrc' or restart your shell to use conda."
fi

# Ensure Conda is active in this shell
export PATH="$HOME/miniconda/bin:$PATH"
eval "$($HOME/miniconda/bin/conda shell.bash hook)"

# --- STEP 4: Create minimal conda environment (optional) ---
echo "[INFO] Skipping large conda environment due to disk space constraints."
echo "[INFO] You can manually create a minimal environment later with:"
echo "       conda env create -f environment-minimal.yml"
echo "[INFO] Or use Docker for full tool support."

# --- STEP 5: Install Nextflow ---
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

echo ""
echo "[INFO] Setup complete! Lightweight demo environment ready."
echo ""
echo "[INFO] To run the SMALL pipeline (quick demo):"
echo "1. Run the simple version check pipeline: nextflow run main.nf --outdir results"
echo "2. This checks tool versions and completes in seconds"
echo ""
echo "[INFO] To run the BIG pipeline (full bioinformatics):"
echo "1. Go to the docker directory: cd ../nextflow-docker-not-supported"
echo "2. Build the image: docker build -t rnaseq-tools:latest ."
echo "3. Run with docker: nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results"
echo ""
echo "[INFO] To use conda in new terminal sessions, run:"
echo "source ~/.bashrc"
echo ""
echo "Or restart your terminal/shell session."
