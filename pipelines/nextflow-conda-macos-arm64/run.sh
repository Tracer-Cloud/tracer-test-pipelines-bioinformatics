#!/bin/bash

# Create logs directory if it doesn't exist
mkdir -p logs

# Set appropriate permissions for logs directory
# 755 = Owner can read/write/execute, Group and Others can read/execute
chmod 755 logs

# Set environment variable for log file location
export NXF_LOG_FILE="logs/nextflow.log"

# Run Nextflow with explicit log file location
nextflow -log logs/nextflow.log run main.nf "$@"