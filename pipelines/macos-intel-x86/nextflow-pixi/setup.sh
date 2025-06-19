#!/bin/bash

# Setup script for Nextflow Pixi pipeline
# This script prepares the environment and validates the setup

set -euo pipefail

echo "=== Nextflow Pixi Pipeline Setup ==="

# Check if Pixi is installed
if ! command -v pixi &> /dev/null; then
    echo "❌ Pixi is not installed!"
    echo "   Install with: curl -fsSL https://pixi.sh/install.sh | bash"
    exit 1
fi

echo "✅ Pixi is available: $(pixi --version)"

# Check if we're in the right directory
if [ ! -f "pixi.toml" ]; then
    echo "❌ pixi.toml not found. Make sure you're in the nextflow-pixi directory"
    exit 1
fi

echo "✅ Found pixi.toml configuration"

# Install dependencies
echo ""
echo "Installing dependencies with Pixi..."
pixi install

# Verify installation
echo ""
echo "Verifying installation..."
pixi run check-env

# Create necessary directories
echo ""
echo "Creating directories..."
mkdir -p results logs
chmod 755 results logs

# Validate test data
echo ""
echo "Validating test data..."
if [ -d "test_data" ] && [ -f "test_data/sample1.fasta" ] && [ -f "test_data/sample2.fasta" ]; then
    echo "✅ Test data is available"
    echo "   - $(wc -l test_data/sample1.fasta | awk '{print $1}') lines in sample1.fasta"
    echo "   - $(wc -l test_data/sample2.fasta | awk '{print $1}') lines in sample2.fasta"
else
    echo "⚠️  Test data not found, creating sample data..."
    mkdir -p test_data
    
    # Create sample FASTA files
    cat > test_data/sample1.fasta << 'EOF'
>sequence1
ATCGATCGATCGATCG
>sequence2
GCTAGCTAGCTAGCTA
>sequence3
TTTTAAAACCCCGGGG
EOF

    cat > test_data/sample2.fasta << 'EOF'
>seq_A
AAAAAAAAA
>seq_B
TTTTTTTTT
>seq_C
CCCCCCCCC
>seq_D
GGGGGGGGG
>seq_E
ATCGATCGAT
EOF
    
    echo "✅ Created sample test data"
fi

echo ""
echo "=== Setup Complete! ==="
echo ""
echo "Next steps:"
echo "1. Run the pipeline: pixi run pipeline"
echo "2. Check results in: results/"
echo "3. View logs in: logs/"
echo ""
echo "Available tasks:"
echo "  pixi run pipeline      # Run main pipeline"
echo "  pixi run test         # Run with test data"
echo "  pixi run clean        # Clean up results"
echo "  pixi run check-env    # Verify environment"
echo "  pixi run dev          # Full development workflow"
