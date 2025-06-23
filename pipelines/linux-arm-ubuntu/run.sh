#!/bin/bash

echo "[INFO] Launching Pixi Shell + Tracer Init"
pixi run tracer init

echo "[INFO] Running Nextflow Pipeline..."
nextflow run main.nf -c nextflow.config
