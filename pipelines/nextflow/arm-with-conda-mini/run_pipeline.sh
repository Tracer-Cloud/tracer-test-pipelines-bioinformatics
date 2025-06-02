#!/bin/bash

# Script to run Nextflow pipeline with organized logging
# Usage: ./run_pipeline.sh [input_pattern]

set -e  # Exit on any error

# Default parameters
INPUT_PATTERN="${1:-test_data/*.fasta}"
LOGS_DIR="logs"
RESULTS_DIR="results"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Nextflow Pipeline Runner ===${NC}"
echo "Input pattern: $INPUT_PATTERN"
echo "Results directory: $RESULTS_DIR"
echo "Logs directory: $LOGS_DIR"
echo ""

# Create directories if they don't exist
echo -e "${YELLOW}Creating directories...${NC}"
mkdir -p "$LOGS_DIR"
mkdir -p "$RESULTS_DIR"

# Check if Java is available
if ! command -v java &> /dev/null; then
    echo -e "${YELLOW}Java not found in PATH. Adding OpenJDK to PATH...${NC}"
    export PATH="/opt/homebrew/opt/openjdk@11/bin:$PATH"
fi

# Check Java version
echo -e "${YELLOW}Java version:${NC}"
java -version

# Set Nextflow log file location
export NXF_LOG_FILE="$LOGS_DIR/.nextflow.log"

# Run the pipeline
echo -e "${GREEN}Running Nextflow pipeline...${NC}"
echo "Command: nextflow run main.nf --input \"$INPUT_PATTERN\" --outdir \"$RESULTS_DIR\" --logsdir \"$LOGS_DIR\""
echo ""

nextflow run main.nf \
    --input "$INPUT_PATTERN" \
    --outdir "$RESULTS_DIR" \
    --logsdir "$LOGS_DIR"

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}=== Pipeline completed successfully! ===${NC}"
    echo ""
    echo -e "${YELLOW}Generated files:${NC}"
    echo "Results: $RESULTS_DIR/"
    echo "Logs: $LOGS_DIR/"
    echo ""
    echo -e "${YELLOW}Log files:${NC}"
    ls -la "$LOGS_DIR/"
    echo ""
    echo -e "${YELLOW}Result files:${NC}"
    find "$RESULTS_DIR" -type f 2>/dev/null || echo "No result files found"
else
    echo ""
    echo -e "${RED}=== Pipeline failed! ===${NC}"
    echo "Check the log files in $LOGS_DIR/ for details"
    exit 1
fi
