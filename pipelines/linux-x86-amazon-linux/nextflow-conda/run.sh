#!/bin/bash
set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}[INFO]${NC} Setting up Nextflow pipeline for Linux x86 Amazon Linux (Conda)..."

# Check available memory and set Java options
AVAILABLE_MEMORY=$(free -m | awk 'NR==2{printf "%.0f", $7}')
echo -e "${BLUE}[INFO]${NC} Available memory: ${AVAILABLE_MEMORY}MB"

if [ "$AVAILABLE_MEMORY" -lt 512 ]; then
    echo -e "${YELLOW}[WARNING]${NC} Low memory system detected. Using minimal memory settings."
    export NXF_OPTS="-Xmx256m"
else
    export NXF_OPTS="-Xmx512m"
fi

# Install Miniconda if not present
if [ ! -d "$HOME/miniconda" ]; then
    echo -e "${BLUE}[INFO]${NC} Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"
    rm miniconda.sh
    
    # Add conda to PATH
    export PATH="$HOME/miniconda/bin:$PATH"
    
    # Initialize conda
    "$HOME/miniconda/bin/conda" init bash
    
    echo -e "${GREEN}[SUCCESS]${NC} Miniconda installed successfully."
else
    echo -e "${BLUE}[INFO]${NC} Miniconda already installed."
    export PATH="$HOME/miniconda/bin:$PATH"
fi

# Ensure conda is properly initialized
if ! conda info --base &> /dev/null; then
    echo -e "${BLUE}[INFO]${NC} Initializing conda..."
    conda init bash
    source ~/.bashrc 2>/dev/null || true
fi

# Check if conda environment exists and setup
if conda env list | grep -q "^rnaseq-minimal[[:space:]]"; then
    echo -e "${BLUE}[INFO]${NC} Conda environment 'rnaseq-minimal' already exists. Updating from environment-minimal.yml..."
    conda env update -f environment-minimal.yml
else
    echo -e "${BLUE}[INFO]${NC} Creating conda environment from environment-minimal.yml..."
    conda env create -f environment-minimal.yml
    echo -e "${GREEN}[SUCCESS]${NC} Environment created successfully!"
fi

# Create necessary directories
mkdir -p test_data logs results test_results

# Create test data if needed
if [ ! -f "test_data/sample1.fasta" ]; then
    echo -e "${BLUE}[INFO]${NC} Creating sample test data..."
    cat > test_data/sample1.fasta << 'EOF'
>seq1
ATCGATCGATCG
>seq2
GCTAGCTAGCTA
EOF
fi

# Run pipeline using conda run
echo -e "${BLUE}[INFO]${NC} Running Nextflow pipeline (version check)..."
echo -e "${BLUE}[INFO]${NC} Using NXF_OPTS: $NXF_OPTS"

if conda run -n rnaseq-minimal nextflow -log logs/nextflow.log run main.nf --outdir test_results; then
    echo -e "${GREEN}[SUCCESS]${NC} Pipeline completed successfully!"
    echo -e "${BLUE}[INFO]${NC} Results available in: test_results/"
    echo -e "${BLUE}[INFO]${NC} Logs available in: logs/nextflow.log"
    
    # Show the version check results
    if [ -f "test_results/tool_versions.txt" ]; then
        echo -e "${BLUE}[INFO]${NC} Version check results:"
        cat test_results/tool_versions.txt
    fi
else
    echo -e "${RED}[ERROR]${NC} Pipeline failed. Check logs/nextflow.log for details"
    exit 1
fi 