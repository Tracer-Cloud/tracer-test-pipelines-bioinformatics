#!/bin/bash

echo "[INFO] Launching Pixi Shell + Tracer Init"

# Check for Pixi, install if missing
if ! command -v pixi &> /dev/null; then
    echo "[ERROR] Pixi not found. Installing Pixi..."
    curl -sSf https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
fi

# Run Pixi init (comment this out if not needed)
pixi run tracer init

# Check for Nextflow, install if missing
if ! command -v nextflow &> /dev/null; then
    echo "[INFO] Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    mv nextflow /usr/local/bin/ 2>/dev/null || sudo mv nextflow /usr/local/bin/
fi

# Confirm Nextflow is available
if ! command -v nextflow &> /dev/null; then
    echo "[ERROR] Nextflow installation failed or is not in PATH"
    exit 1
fi

echo "[INFO] Running Nextflow Pipeline..."
nextflow run main.nf -c nextflow.config
