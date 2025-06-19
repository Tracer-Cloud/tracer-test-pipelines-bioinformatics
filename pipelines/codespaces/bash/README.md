# Codespaces Pixi Pipeline Runner

This directory contains a bash script that sets up and runs a Nextflow pipeline using Pixi in GitHub Codespaces.

## 🚀 Quick Start

```bash
# Navigate to this directory
cd pipelines/codespaces/bash

# Run the setup and pipeline script
./run_pixi_pipeline.sh
```

## 📁 What it does

The script will:
1. Install Pixi if not already installed
2. Create a Linux-compatible version of the macOS ARM64 Pixi pipeline
3. Set up the environment and dependencies
4. Run the Nextflow pipeline
5. Display results

## 🔧 Features

- **Cross-platform**: Adapts the macOS Pixi example for Linux (Codespaces)
- **Automatic setup**: Installs all dependencies automatically
- **Test data**: Includes sample FASTA files for testing
- **Clean output**: Shows pipeline results and logs

## 🧬 Pipeline Description

The pipeline processes FASTA files and performs:
1. **Version Check**: Displays versions of FastQC and STAR
2. **FASTA Statistics**: Calculates sequence counts, lengths, and statistics
3. **Sequence Counting**: Counts sequences in each FASTA file

## 📊 Expected Output

After running, you'll find:
- `results/fasta_stats/`: Detailed statistics for each FASTA file
- `results/counts/`: Sequence counts for each file
- `logs/`: Pipeline execution logs and reports