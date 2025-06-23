#!/bin/bash

# Step 0: Ensure Nextflow is installed
if ! command -v nextflow &> /dev/null
then
    echo "[INFO] Nextflow not found, installing..."
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

# Step 1: Run the pipeline
echo "[INFO] Running Nextflow pipeline..."
nextflow run main.nf -c nextflow.config

