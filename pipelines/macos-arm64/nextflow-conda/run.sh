#!/usr/bin/env bash
set -euo pipefail

echo "Running Nextflow pipeline..."

# Install Homebrew if not present
if ! command -v brew &> /dev/null; then
    echo "Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    
    # Add Homebrew to PATH
    if [[ $(uname -m) == "arm64" ]]; then
        eval "$(/opt/homebrew/bin/brew shellenv)"
    else
        eval "$(/usr/local/bin/brew shellenv)"
    fi
fi

# Install conda via Homebrew if not present
if ! command -v conda &> /dev/null; then
    echo "Installing Miniforge via Homebrew..."
    brew install miniforge
    conda init "$(basename "$SHELL")"
    source ~/.zshrc 2>/dev/null || source ~/.bashrc 2>/dev/null
fi

# Ensure conda is properly initialized for current shell
if ! conda info --base &> /dev/null; then
    echo "Initializing conda for current shell..."
    conda init "$(basename "$SHELL")"
    # Source the shell configuration
    if [[ "$SHELL" == *"zsh"* ]]; then
        source ~/.zshrc
    elif [[ "$SHELL" == *"bash"* ]]; then
        source ~/.bashrc
    fi
fi

# Setup conda environment if not present
if ! conda env list | grep -q "nextflow-minimal"; then
    echo "Creating conda environment 'nextflow-minimal'..."
    conda env create -f environment.yml
fi

echo "Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow-minimal

# Create directories if they don't exist
mkdir -p test_data logs results test_results

# Create sample test data if it doesn't exist
if [ ! -f "test_data/sample1.fasta" ]; then
    echo "Creating sample test data..."
    cat > test_data/sample1.fasta << 'EOF'
>seq1
ATCGATCGATCG
>seq2
GCTAGCTAGCTA
EOF
fi

# Run the Nextflow pipeline
nextflow run main.nf --input test_data/*.fasta --iterations 20 --outdir test_results