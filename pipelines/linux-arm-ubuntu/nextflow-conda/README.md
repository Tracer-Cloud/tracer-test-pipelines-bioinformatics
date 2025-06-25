# Nextflow Pipeline on Linux ARM Ubuntu - Conda Version

This directory contains a **Conda-based Nextflow pipeline** for ARM Ubuntu. The setup creates a dedicated conda environment with all necessary dependencies.

## Quick Start (Recommended: Conda)

### Prerequisites

- Ubuntu ARM64 system
- Internet connection for downloading dependencies

### Initial Setup

```bash
# Run the setup script to install conda and create the environment (if needed)
./run.sh
```

The setup script will:

1. Install Miniconda for ARM64 (if not present)
2. Create the conda environment named `linux-arm-ubuntu-minimal` with all dependencies (if not already present)
3. Activate the environment and run the pipeline

> **Note:** If the environment `linux-arm-ubuntu-minimal` already exists, the script will use it as is and will not recreate or update it.

### Using the Pipeline

After running `run.sh`, the pipeline will execute automatically in the correct environment.

#### Manual Usage

If you want to activate the environment and run commands manually:

```bash
# Activate the conda environment
conda activate linux-arm-ubuntu-minimal

# Run the basic pipeline
nextflow run main.nf --outdir results

# Deactivate when done
conda deactivate
```
