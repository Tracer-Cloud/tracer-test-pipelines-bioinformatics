#!/bin/bash

# Test script to verify the Codespaces Pixi setup works
# This is a lightweight test that doesn't run the full pipeline

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${BLUE}[TEST]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

print_error() {
    echo -e "${RED}[FAIL]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

print_status "Testing Codespaces Pixi Pipeline Setup"

# Test 1: Check if script exists and is executable
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_SCRIPT="$SCRIPT_DIR/run_pixi_pipeline.sh"

if [ -f "$MAIN_SCRIPT" ] && [ -x "$MAIN_SCRIPT" ]; then
    print_success "Main script exists and is executable"
else
    print_error "Main script not found or not executable: $MAIN_SCRIPT"
    exit 1
fi

# Test 2: Check if we can install Pixi (dry run)
print_status "Checking Pixi installation availability..."
if curl -fsSL https://pixi.sh/install.sh | head -10 | grep -q "pixi"; then
    print_success "Pixi installer is accessible"
else
    print_error "Cannot access Pixi installer"
    exit 1
fi

# Test 3: Check if we can access Nextflow installer
print_status "Checking Nextflow installation availability..."
if curl -s --max-time 10 https://get.nextflow.io | head -10 | grep -q "nextflow" 2>/dev/null; then
    print_success "Nextflow installer is accessible"
else
    print_status "Nextflow installer check skipped (network timeout or unavailable)"
fi

# Test 4: Check if macOS pipeline source exists
WORKSPACE_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
MACOS_PIPELINE="$WORKSPACE_ROOT/pipelines/macos-arm64/nextflow-pixi"

if [ -d "$MACOS_PIPELINE" ]; then
    print_success "macOS ARM64 pipeline source found"
    
    # Check for required files
    required_files=("main.nf" "nextflow.config" "run.sh" "pixi.toml")
    for file in "${required_files[@]}"; do
        if [ -f "$MACOS_PIPELINE/$file" ]; then
            print_success "Required file found: $file"
        else
            print_error "Required file missing: $file"
            exit 1
        fi
    done
else
    print_error "macOS ARM64 pipeline source not found: $MACOS_PIPELINE"
    exit 1
fi

# Test 5: Verify script syntax
print_status "Checking script syntax..."
if bash -n "$MAIN_SCRIPT"; then
    print_success "Script syntax is valid"
else
    print_error "Script has syntax errors"
    exit 1
fi

print_success "All tests passed! The Codespaces Pixi setup should work correctly."
print_status "To run the full pipeline, execute: ./run_pixi_pipeline.sh"
