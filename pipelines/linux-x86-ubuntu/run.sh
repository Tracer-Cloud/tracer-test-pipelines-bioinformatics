#!/bin/bash
set -e

echo "[INFO] Starting environment setup for linux-x86-ubuntu..."

# Install Java if not present
if ! command -v java &> /dev/null; then
    echo "[INFO] Installing Java..."
    apt-get update && apt-get install -y openjdk-17-jdk
else
    echo "[INFO] Java is already installed."
fi

# Install Miniconda if not present
if [ ! -d "$HOME/miniconda" ]; then
    echo "[INFO] Installing Miniconda for x86_64..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    rm miniconda.sh
else
    echo "[INFO] Miniconda already installed."
fi

export PATH="$HOME/miniconda/bin:$PATH"

# Initialize conda if needed
if ! grep -q "conda initialize" ~/.bashrc; then
    $HOME/miniconda/bin/conda init bash
    source ~/.bashrc
fi

# Install Nextflow if not present
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

mkdir -p test_data results

echo "[INFO] Setup complete."
echo "Java version:"
java -version
echo "Conda version:"
conda --version
echo "Nextflow version:"
nextflow -version
