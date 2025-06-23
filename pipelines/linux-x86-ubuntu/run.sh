#!/bin/bash
set -e

echo "[INFO] Starting environment setup for linux-x86-ubuntu..."

# Step 1: Install Java
if ! command -v java &>/dev/null; then
    echo "[INFO] Installing OpenJDK 17..."
    apt-get update
    apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java is already installed."
fi

# Step 2: Install Miniconda (fresh, safe)
CONDA_DIR="/root/miniconda"
if [ ! -x "$CONDA_DIR/bin/conda" ]; then
    echo "[INFO] Installing Miniconda for x86_64..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    rm -rf "$CONDA_DIR"
    bash miniconda.sh -b -p "$CONDA_DIR"
else
    echo "[INFO] Miniconda already installed."
fi

# Make sure conda is available in this shell session
export PATH="$CONDA_DIR/bin:$PATH"

# Step 3: Install Nextflow
if ! command -v nextflow &>/dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

# Final test output
echo -e "\n[✅] Setup complete."
echo "🔍 Java:"
java -version
echo "🔍 Conda:"
conda --version || echo "❌ Conda still not found in PATH"
echo "🔍 Nextflow:"
nextflow -version || echo "❌ Nextflow not found"
