#!/bin/bash

# Step 0: Ensure Java is installed
if ! command -v java &> /dev/null; then
    echo "[INFO] Java not found. Installing OpenJDK 17..."
    apt-get update -y
    apt-get install -y openjdk-17-jre
else
    echo "[INFO] Java already installed."
fi

# Step 1: Ensure Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Nextflow not found, installing..."
    curl -s https://get.nextflow.io | bash
    mv nextflow /usr/local/bin/
else
    echo "[INFO] Nextflow already installed."
fi

# Step 2: Run the pipeline
echo "[INFO] Running Nextflow pipeline..."
nextflow run main.nf -c nextflow.config
