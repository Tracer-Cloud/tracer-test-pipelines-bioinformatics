#!/bin/bash

echo "[INFO] Launching Pixi Shell + Tracer Init"

if ! command -v pixi &> /dev/null; then
    echo "[ERROR] Pixi not found. Please install Pixi first:"
    echo "  curl -sSf https://pixi.sh/install.sh | bash"
    exit 1
fi

pixi run tracer init  # if needed, or comment/remove if not required

echo "[INFO] Running Nextflow Pipeline..."
nextflow run main.nf -c nextflow.config
