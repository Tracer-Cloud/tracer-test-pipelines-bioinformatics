#!/usr/bin/env bash
set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}[INFO]${NC} Setting up Nextflow pipeline for Linux x86 Amazon Linux..."

# Check available memory
AVAILABLE_MEMORY=$(free -m | awk 'NR==2{printf "%.0f", $7}')
echo -e "${BLUE}[INFO]${NC} Available memory: ${AVAILABLE_MEMORY}MB"

if [ "$AVAILABLE_MEMORY" -lt 512 ]; then
    echo -e "${YELLOW}[WARNING]${NC} Low memory system detected. Using minimal memory settings."
    export NXF_OPTS="-Xmx256m"
else
    export NXF_OPTS="-Xmx512m"
fi

# Install Pixi if not present
if ! command -v pixi &> /dev/null; then
    echo -e "${BLUE}[INFO]${NC} Installing Pixi..."
    curl -fsSL https://pixi.sh/install.sh | bash
    export PATH="$HOME/.pixi/bin:$PATH"
fi

# Source shell profile if available
if [ -n "$ZSH_VERSION" ]; then
    SHELL_PROFILE="$HOME/.zshrc"
elif [ -n "$BASH_VERSION" ]; then
    SHELL_PROFILE="$HOME/.bashrc"
else
    SHELL_PROFILE="$HOME/.profile"
fi

if [ -f "$SHELL_PROFILE" ]; then
    echo -e "${BLUE}[INFO]${NC} Sourcing $SHELL_PROFILE..."
    source "$SHELL_PROFILE"
fi

# Install dependencies
echo -e "${BLUE}[INFO]${NC} Installing Pixi dependencies..."
pixi install

# Create necessary directories
mkdir -p logs results test_results

# Run environment check
echo -e "${BLUE}[INFO]${NC} Running environment check..."
if pixi run check-env; then
    echo -e "${GREEN}[SUCCESS]${NC} Environment check passed!"
else
    echo -e "${YELLOW}[WARNING]${NC} Environment check had issues, but continuing..."
fi

# Run pipeline with memory limits
echo -e "${BLUE}[INFO]${NC} Running Nextflow pipeline..."
echo -e "${BLUE}[INFO]${NC} Using NXF_OPTS: $NXF_OPTS"

if pixi run pipeline; then
    echo -e "${GREEN}[SUCCESS]${NC} Pipeline completed successfully!"
    echo -e "${BLUE}[INFO]${NC} Results available in: results/"
    echo -e "${BLUE}[INFO]${NC} Logs available in: logs/nextflow.log"
else
    echo -e "${RED}[ERROR]${NC} Pipeline failed. Check logs/nextflow.log for details"
    exit 1
fi
