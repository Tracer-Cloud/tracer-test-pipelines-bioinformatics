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

## Files

- `main.nf` - Nextflow pipeline definition (shared from `../../shared/nextflow/`)
- `nextflow.config` - Nextflow configuration (shared from `../../shared/nextflow/`)
- `custom.config` - Custom configuration for nf-core pipelines (shared from `../../shared/nextflow/`)
- `test_data/` - Sample test data (shared from `../../shared/nextflow/test_data/`)
- `pixi.toml` - Pixi configuration with dependencies and tasks
- `pixi.lock` - Pixi lock file for reproducible builds

## Notes

- Test data is referenced from the shared folder to maintain consistency across all platforms
- The pipeline uses explicit version specifiers to avoid Pixi warnings
- The nf-core RNA-seq pipeline uses a specific revision (3.14.0) to avoid version conflicts
