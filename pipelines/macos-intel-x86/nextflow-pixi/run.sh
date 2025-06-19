#!/bin/bash

# Nextflow pipeline runner for Pixi environment
# This script assumes it's being run within a Pixi environment

set -euo pipefail

echo "=== Nextflow Pipeline with Pixi Environment ==="
echo "Starting pipeline execution..."

# Create logs directory if it doesn't exist
mkdir -p logs

# Set appropriate permissions for logs directory
chmod 755 logs

# Set environment variable for log file location
export NXF_LOG_FILE="logs/nextflow.log"

# Verify we're in a Pixi environment
if command -v pixi &> /dev/null; then
    echo "✅ Pixi is available"
else
    echo "⚠️  Pixi not found - make sure to run this within a Pixi environment"
    echo "   Use: pixi run pipeline"
fi

# Check if required tools are available
echo "Checking tool availability..."
if command -v fastqc &> /dev/null; then
    echo "✅ FastQC available: $(fastqc --version 2>&1 | head -1)"
else
    echo "❌ FastQC not found"
    exit 1
fi

if command -v STAR &> /dev/null; then
    echo "✅ STAR available: $(STAR --version 2>&1 | head -1)"
else
    echo "❌ STAR not found"
    exit 1
fi

if command -v nextflow &> /dev/null; then
    echo "✅ Nextflow available: $(nextflow -v 2>&1 | head -1)"
else
    echo "❌ Nextflow not found"
    echo "   Install with: pixi add nextflow"
    exit 1
fi

echo ""
echo "=== Running Nextflow Pipeline ==="

# Run Nextflow with explicit log file location
# Pass all command line arguments to nextflow
nextflow -log logs/nextflow.log run main.nf "$@"

echo ""
echo "=== Pipeline Completed ==="
echo "Check results in: results/"
echo "Check logs in: logs/"
