# Nextflow Pipeline on Linux ARM Ubuntu - Pixi Version

This directory contains a **Pixi-based Nextflow pipeline** for ARM Ubuntu. The setup and usage are similar to the shared Nextflow pipeline, but adapted for ARM Ubuntu.

## Quick Start (Recommended: Pixi)

### Prerequisites

- [Pixi](https://pixi.sh) installed on your system
- ARM Ubuntu system

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
