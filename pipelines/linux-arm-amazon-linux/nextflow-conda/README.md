# Nextflow Pipeline on Linux ARM Amazon Linux - Conda Version

This directory contains a **Conda-based Nextflow pipeline** for ARM Amazon Linux. The setup and usage are similar to the x86 Ubuntu reference, but adapted for ARM and Amazon Linux.

## Quick Start - Conda Pipeline

### Prerequisites

- ARM Amazon Linux system
- sudo access for package installation

### One-Command Setup

1. **Run the setup script:**

   ```bash
   ./run.sh
   ```

2. **Activate conda in your current shell session:**

   ```bash
   source ~/.bashrc
   ```

3. **Verify conda is available:**
   ```bash
   conda --version
   ```

## Files

- `main.nf` - Nextflow pipeline definition
- `nextflow.config` - Nextflow configuration
- `custom.config` - Custom configuration for nf-core pipelines
- `test_data/` - Sample test data for the pipeline
- `run.sh` - Setup script for ARM Amazon Linux

## Usage

### Run the pipeline

```bash
nextflow run main.nf --outdir results
```

### Run nf-core RNA-seq pipeline

```bash
nextflow run nf-core/rnaseq -c custom.config -profile test --outdir results -resume
```

## Troubleshooting

If you encounter issues with conda or nextflow, ensure you have sourced your shell profile or restarted your terminal after running the setup script.

---

For more details, see the x86 Ubuntu README in `../linux-x86-ubuntu/nextflow-conda/README.md`.
