# Nextflow Pipeline on Linux ARM Ubuntu - Conda Version

This directory contains a **Conda-based Nextflow pipeline** for ARM Ubuntu. The setup and usage are similar to the shared Nextflow pipeline, but adapted for ARM Ubuntu.

## Quick Start - Conda Pipeline

### Prerequisites

- ARM Ubuntu system
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

- `main.nf` - Nextflow pipeline definition (shared from `../../shared/nextflow/`)
- `nextflow.config` - Nextflow configuration (shared from `../../shared/nextflow/`)
- `custom.config` - Custom configuration for nf-core pipelines (shared from `../../shared/nextflow/`)
- `test_data/` - Sample test data (shared from `../../shared/nextflow/test_data/`)
- `run.sh` - Setup script for ARM Ubuntu

## Usage

### Run the pipeline

```bash
nextflow run main.nf --outdir results
```

### Run nf-core RNA-seq pipeline

```bash
nextflow run nf-core/rnaseq -r 3.14.0 -c custom.config -profile test --outdir results -resume
```

## Troubleshooting

If you encounter issues with conda or nextflow, ensure you have sourced your shell profile or restarted your terminal after running the setup script.

## Notes

- Test data is referenced from the shared folder to maintain consistency across all platforms
- The nf-core RNA-seq pipeline uses a specific revision (3.14.0) to avoid version conflicts

---

For more details, see the shared README in `../../shared/nextflow/README.md`.
