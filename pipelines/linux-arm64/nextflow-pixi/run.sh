#!/bin/bash

# Run script for Nextflow Pixi pipeline (Linux ARM64)
# This script runs the pipeline using Pixi-managed dependencies

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Default parameters
INPUT_PATTERN="test_data/*.fasta"
OUTPUT_DIR="results"
PROFILE="standard"
RESUME=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_PATTERN="$2"
            shift 2
            ;;
        --outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --resume)
            RESUME="-resume"
            shift
            ;;
        --help)
            echo "Usage: $0 [--input PATTERN] [--outdir DIR] [--profile PROFILE] [--resume] [--help]"
            echo "  --input    Input FASTA files pattern (default: test_data/*.fasta)"
            echo "  --outdir   Output directory (default: results)"
            echo "  --profile  Execution profile (default: standard)"
            echo "  --resume   Resume previous run"
            echo "  --help     Show this help"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

print_status "=== Nextflow Pixi Pipeline Runner (Linux ARM64) ==="
print_status "Input: $INPUT_PATTERN"
print_status "Output: $OUTPUT_DIR"
print_status "Profile: $PROFILE"

# Check if we're in a Pixi environment
if [ -z "${PIXI_PROJECT_ROOT:-}" ]; then
    print_error "Not in a Pixi environment. Please run 'pixi shell' first or use 'pixi run pipeline'"
    exit 1
fi

# Create directories
mkdir -p "$OUTPUT_DIR" logs

# Check input files
if ! ls $INPUT_PATTERN >/dev/null 2>&1; then
    print_error "No input files found matching: $INPUT_PATTERN"
    exit 1
fi

print_success "Found input files: $(ls $INPUT_PATTERN | wc -l)"

# Run the pipeline
print_status "Running Nextflow pipeline..."
nextflow run main.nf \
    --input "$INPUT_PATTERN" \
    --outdir "$OUTPUT_DIR" \
    -profile "$PROFILE" \
    $RESUME \
    -with-trace logs/trace.txt \
    -with-report logs/report.html \
    -with-timeline logs/timeline.html

if [ $? -eq 0 ]; then
    print_success "Pipeline completed successfully!"
    print_status "Results available in: $OUTPUT_DIR"
    print_status "Reports available in: logs/"
else
    print_error "Pipeline failed!"
    exit 1
fi
