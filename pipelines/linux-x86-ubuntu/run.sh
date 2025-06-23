#!/bin/bash
set -e

echo "[INFO] Setting up environment..."

# Java
if ! command -v java &>/dev/null; then
  echo "[INFO] Installing OpenJDK..."
  sudo apt-get update
  sudo apt-get install -y default-jdk
fi

# Conda
if ! command -v conda &>/dev/null; then
  echo "[INFO] Installing Miniconda for x86_64..."
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p "$HOME/miniconda"
  export PATH="$HOME/miniconda/bin:$PATH"
  echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
  source ~/.bashrc
fi

# Nextflow
if ! command -v nextflow &>/dev/null; then
  echo "[INFO] Installing Nextflow..."
  curl -s https://get.nextflow.io | bash
  chmod +x nextflow
  sudo mv nextflow /usr/local/bin/
fi

echo "[✅] Environment ready."
