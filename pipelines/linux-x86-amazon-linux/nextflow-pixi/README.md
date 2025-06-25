# Nextflow Pipeline on Linux x86 Amazon Linux - Pixi Version

This directory contains a **Pixi-based Nextflow pipeline** for x86 Amazon Linux. The setup and usage are similar to the x86 Ubuntu reference, but adapted for Amazon Linux.

## Quick Start (Recommended: Pixi)

### Prerequisites

- [Pixi](https://pixi.sh) installed on your system
- x86 Amazon Linux system

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
