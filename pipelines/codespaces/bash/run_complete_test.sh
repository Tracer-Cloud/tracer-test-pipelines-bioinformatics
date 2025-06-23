#!/bin/bash
# Complete test runner - executes all steps in order
set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log "=== Bioinformatics Pipeline Complete Test ==="

# Check if we're in the right directory
if [ ! -f "run.sh" ]; then
    echo "Error: Please run this script from the pipelines/codespaces/bash directory"
    exit 1
fi

# Step 1: Make all scripts executable
log "Step 1: Making scripts executable..."
chmod +x *.sh

# Step 2: Setup reference genome
log "Step 2: Setting up reference genome..."
if [ ! -d "reference" ] || [ ! -f "reference/chr19.fa" ]; then
    ./setup_reference.sh
else
    log "Reference genome already exists, skipping setup"
fi

# Step 3: Generate test data (if needed)
log "Step 3: Checking test data..."
if [ ! -f "test_data/test_1.fastq" ] || [ ! -f "test_data/test_2.fastq" ]; then
    log "Generating test data..."
    ./generate_test_data.sh
else
    log "Test data already exists, skipping generation"
fi

# Step 4: Run the pipeline with existing test files
log "Step 4: Running pipeline with existing test files..."
if [ -f "small_test_fastq_1.fastq" ] && [ -f "small_test_fastq_2.fastq" ]; then
    ./run.sh small_test_fastq_1.fastq small_test_fastq_2.fastq output/test_output.vcf
else
    log "Using generated test files..."
    ./run.sh test_data/test_1.fastq test_data/test_2.fastq output/test_output.vcf
fi

log "=== Test completed successfully! ==="
log "Check the output directory for results:"
ls -la output/ 