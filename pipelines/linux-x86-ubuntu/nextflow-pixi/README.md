# Nextflow Pipeline on Linux x86 using Conda

# Linux x86-64 Ubuntu Nextflow Pipeline - Pixi Version

This directory contains a **Pixi-based Nextflow pipeline** for modern, fast package management. Pixi provides faster dependency resolution and better reproducibility than traditional conda.

## Quick Start (Recommended: Pixi)

### Prerequisites

- [Pixi](https://pixi.sh) installed on your system
- Linux x86-64 system

### Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Run the Pipeline with Pixi

```bash
# Install dependencies
pixi install

# Check environment
pixi run check-env

# Run the pipeline
pixi run pipeline

# Run with custom output directory
pixi run rnaseq

# Clean up
pixi run clean
```
