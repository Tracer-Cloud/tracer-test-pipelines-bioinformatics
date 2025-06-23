#!/bin/bash

# Step 0: Ensure Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Nextflow not found. Installing..."

    # Download and make executable
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow

    # Move to PATH and rehash
    mv nextflow /usr/local/bin/nextflow
    export PATH=$PATH:/usr/local/bin

    echo "[INFO] Nextflow installed to /usr/local/bin"
else
    echo "[INFO] Nextflow already installed."
fi

# Step 1: Run the pipeline
echo "[INFO] Running Nextflow pipeline..."
nextflow run main.nf -c nextflow.config
