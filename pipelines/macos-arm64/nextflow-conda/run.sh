#!/usr/bin/env bash
set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}[INFO]${NC} Setting up Nextflow pipeline for macOS (Conda only)..."

# Install Homebrew if not present
if ! command -v brew &> /dev/null; then
    echo -e "${BLUE}[INFO]${NC} Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    
    # Add Homebrew to PATH
    if [[ $(uname -m) == "arm64" ]]; then
        eval "$(/opt/homebrew/bin/brew shellenv)"
    else
        eval "$(/usr/local/bin/brew shellenv)"
    fi
fi

# Install conda via Homebrew
if ! command -v conda &> /dev/null; then
    echo -e "${BLUE}[INFO]${NC} Installing Miniforge via Homebrew..."
    brew install miniforge
    conda init "$(basename "$SHELL")"
    source ~/.zshrc 2>/dev/null || source ~/.bashrc 2>/dev/null
fi

# Ensure conda is properly initialized for current shell
if ! conda info --base &> /dev/null; then
    echo -e "${BLUE}[INFO]${NC} Initializing conda for current shell..."
    conda init "$(basename "$SHELL")"
    # Source the shell configuration
    if [[ "$SHELL" == *"zsh"* ]]; then
        source ~/.zshrc
    elif [[ "$SHELL" == *"bash"* ]]; then
        source ~/.bashrc
    fi
fi

# Setup conda environment
echo -e "${BLUE}[INFO]${NC} Setting up conda environment..."
if ! conda env list | grep -q "nextflow-minimal"; then
    echo -e "${BLUE}[INFO]${NC} Creating conda environment 'nextflow-minimal'..."
    conda env create -f environment.yml
else
    echo -e "${BLUE}[INFO]${NC} Updating existing conda environment 'nextflow-minimal'..."
    conda env update -f environment.yml
fi

# Create test data if needed
mkdir -p test_data logs results test_results
if [ ! -f "test_data/sample1.fasta" ]; then
    echo -e "${BLUE}[INFO]${NC} Creating sample test data..."
    cat > test_data/sample1.fasta << 'EOF'
>seq1
ATCGATCGATCG
>seq2
GCTAGCTAGCTA
EOF
fi

# Activate conda environment and run pipeline
echo -e "${BLUE}[INFO]${NC} Activating conda environment and running Nextflow pipeline..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow-minimal

if nextflow -log logs/nextflow.log run main.nf --input test_data/*.fasta --outdir test_results; then
    echo -e "${GREEN}[SUCCESS]${NC} Pipeline completed! Results in: test_results/"
    echo -e "${BLUE}[INFO]${NC} Logs available in: logs/nextflow.log"
else
    echo -e "${RED}[ERROR]${NC} Pipeline failed. Check logs/nextflow.log"
    exit 1
fi
