#!/bin/bash

# Linux x86-64 Ubuntu Pipeline Setup
# This script sets up Conda, Pixi, and Nextflow for bioinformatics pipelines
# Supports both traditional Conda and modern Pixi package management

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_status "Minimal setup for Ubuntu x86_64..."

# --- STEP 1: Install Pixi ---
print_status "Checking Pixi installation..."
if ! command_exists pixi; then
    print_status "Installing Pixi..."
    curl -fsSL https://pixi.sh/install.sh | bash

    # Add Pixi to PATH for current session
    export PATH="$HOME/.pixi/bin:$PATH"

    # Add Pixi to .bashrc if not already present
    if ! grep -q '.pixi/bin' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> "$HOME/.bashrc"
    fi

    # Verify installation
    if command_exists pixi; then
        print_success "Pixi installed successfully"
        pixi --version
    else
        print_error "Failed to install Pixi"
        exit 1
    fi
else
    print_success "Pixi is already installed"
    pixi --version
fi

# --- STEP 2: Install Java ---
if ! command_exists java; then
    print_status "Installing OpenJDK 17..."
    sudo apt-get update && sudo apt-get install -y openjdk-17-jdk
    print_success "Java installed successfully"
else
    print_success "Java already installed"
fi

# --- STEP 3: Install Miniconda (x86_64) ---
if [ ! -d "$HOME/miniconda" ]; then
    print_status "Installing Miniconda (x86_64)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p "$HOME/miniconda"

    # Add conda to PATH in .bashrc if not already present
    if ! grep -q 'miniconda/bin' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> "$HOME/.bashrc"
    fi

    # Initialize conda for bash shell
    "$HOME/miniconda/bin/conda" init bash

    print_success "Miniconda installed successfully"
    print_status "Please run 'source ~/.bashrc' or restart your shell to use conda."
else
    print_success "Miniconda already installed"
fi

# Ensure Conda is active in this shell
export PATH="$HOME/miniconda/bin:$PATH"
eval "$($HOME/miniconda/bin/conda shell.bash hook)"

# --- STEP 4: Install Nextflow ---
if ! command_exists nextflow; then
    print_status "Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    sudo mv nextflow /usr/local/bin/
    print_success "Nextflow installed successfully"
else
    print_success "Nextflow already installed"
fi

# --- Final check ---
print_success "Environment setup complete!"
echo ""
print_status "Installed versions:"
echo "Java version:"
java -version
echo ""
echo "Conda version:"
conda --version
echo ""
echo "Nextflow version:"
nextflow -version
echo ""
echo "Pixi version:"
pixi --version

echo ""
print_status "Usage options:"
echo "1. Use Pixi (recommended): pixi install && pixi run pipeline"
echo "2. Use Conda: source ~/.bashrc && nextflow run main.nf"
echo ""
print_warning "To use conda in new terminal sessions, run: source ~/.bashrc"
print_warning "Or restart your terminal/shell session."
