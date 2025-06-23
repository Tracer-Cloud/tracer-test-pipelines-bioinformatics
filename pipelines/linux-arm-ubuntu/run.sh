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

# Step 2: Prompt for user to run `tracer init`
echo -e "\n🔧 Now run: tracer init"
echo "🛑 The pipeline will stop here so you can configure Tracer."

read -p "👉 Press Enter AFTER you've run 'tracer init' to continue..."

# Step 3: Run the pipeline (user is assumed to have edited config properly)
echo "[INFO] Running Nextflow pipeline..."
nextflow run main.nf -c nextflow.config
