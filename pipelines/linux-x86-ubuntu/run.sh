#!/bin/bash
set -e

echo "[INFO] Starting environment setup for linux-x86-ubuntu..."

# Java
if ! command -v java &>/dev/null; then
  echo "[INFO] Installing Java..."
  apt-get update && apt-get install -y default-jdk
else
  echo "[INFO] Java is already installed."
fi

# Miniconda
if ! command -v conda &>/dev/null; then
  echo "[INFO] Installing Miniconda for x86_64..."
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p $HOME/miniconda
  export PATH="$HOME/miniconda/bin:$PATH"
else
  echo "[INFO] Miniconda already installed."
fi

# Nextflow
if ! command -v nextflow &>/dev/null; then
  echo "[INFO] Installing Nextflow..."
  curl -s https://get.nextflow.io | bash
  mv nextflow /usr/local/bin/
else
  echo "[INFO] Nextflow already installed."
fi

# Run pipeline
echo "[INFO] Running Nextflow pipeline..."
nextflow run main.nf -c nextflow.config
