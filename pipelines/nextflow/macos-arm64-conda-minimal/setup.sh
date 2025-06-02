#!/bin/bash

set -e

echo "🧬 Setting up Nextflow Conda Pipeline on Mac"
echo "=============================================="

# Check if conda/mamba is installed
if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "❌ Conda/Mamba not found. Please install Miniconda or Anaconda first:"
    echo "   brew install miniconda"
    echo "   or visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "📦 Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    echo "✅ Nextflow installed"
else
    echo "✅ Nextflow already installed: $(nextflow -version | head -1)"
fi

# Install mamba if not present (faster than conda)
if ! command -v mamba &> /dev/null; then
    echo "📦 Installing Mamba for faster package resolution..."
    conda install -n base -c conda-forge mamba -y
fi

# Create test data
echo "📄 Creating test FASTA files..."
mkdir -p test_data

# Create sample FASTA file 1
cat > test_data/sample1.fasta << 'EOF'
>sequence1
ATCGATCGATCGATCG
>sequence2
GCTAGCTAGCTAGCTA
>sequence3
TTTTAAAACCCCGGGG
EOF

# Create sample FASTA file 2
cat > test_data/sample2.fasta << 'EOF'
>seq_A
AAAAAAAAAA
>seq_B
TTTTTTTTTT
>seq_C
CCCCCCCCCC
>seq_D
GGGGGGGGGG
>seq_E
ATCGATCGAT
EOF

echo "✅ Test data created in test_data/"

# Create results directory
mkdir -p results

echo ""
echo "🚀 Setup complete! Now you can run the pipeline:"
echo ""
echo "   # Run with test data"
echo "   nextflow run main.nf --input 'test_data/*.fasta'"
echo ""
echo "   # Run with Mac profile (more resources)"
echo "   nextflow run main.nf --input 'test_data/*.fasta' -profile mac"
echo ""
echo "   # Run with custom output directory"
echo "   nextflow run main.nf --input 'test_data/*.fasta' --outdir my_results"
echo ""
echo "📊 Results will be saved in the 'results' directory"
echo "📈 Pipeline reports will be generated automatically"