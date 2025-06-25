# Nextflow Pipeline on Linux ARM Ubuntu - Conda Version

This directory contains a **Conda-based Nextflow pipeline** for ARM Ubuntu. The setup creates a dedicated conda environment with all necessary dependencies.

## Quick Start (Recommended: Conda)

### Prerequisites

- Ubuntu ARM64 system
- Internet connection for downloading dependencies

### Initial Setup

```bash
# Run the setup script to install conda and create environment
chmod +x run.sh
./run.sh
```

The setup script will:

1. Install Docker (if not present)
2. Install OpenJDK 17 (if not present)
3. Install Miniconda for ARM64
4. Create a conda environment named `nextflow-pipeline` with all dependencies
5. Install Nextflow in the conda environment

### Using the Pipeline

#### Option 1: Automatic Pipeline Runner (Recommended)

```bash
# Make the script executable and run
chmod +x run_pipeline.sh
./run_pipeline.sh
```

This script automatically:

- Activates the conda environment
- Verifies Nextflow is available
- Runs the pipeline
- Handles errors gracefully

#### Option 2: Manual Commands

After setup, you need to activate the conda environment before running any commands:

```bash
# Activate the conda environment

conda env create -f environment-minimal.yml
conda activate nextflow-pipeline

# Verify the environment
conda list

# Run the basic pipeline
nextflow run main.nf --outdir results

# Run nf-core RNA-seq pipeline
nextflow run nf-core/rnaseq -r 3.14.0 -c custom.config -profile test --outdir results -resume

# Deactivate when done
conda deactivate
```

### Important Notes

- **Always activate the conda environment** before running Nextflow commands
- The environment is named `nextflow-pipeline`
- If you open a new terminal, you'll need to activate the environment again
- To make conda available in new terminals, run: `source ~/.bashrc`

### Environment Details

The conda environment includes:

- Nextflow
- Java (OpenJDK 17)
- All bioinformatics tools specified in `environment-minimal.yml`

### Troubleshooting

If you encounter issues:

1. **Conda not found**: Run `source ~/.bashrc` or restart your terminal
2. **Environment not found**: Re-run `./run.sh` to recreate the environment
3. **Permission issues**: Ensure Docker is properly configured and your user is in the docker group

### Alternative: Manual Setup

If you prefer to set up manually:

```bash
# Install conda dependencies
conda env create -f environment-minimal.yml

# Activate environment
conda activate nextflow-pipeline

# Install Nextflow
conda install -c bioconda nextflow -y
```

## File Structure

- `main.nf` - Main Nextflow pipeline (links to shared version)
- `nextflow.config` - Nextflow configuration (links to shared version)
- `custom.config` - Custom configuration for this environment
- `environment-minimal.yml` - Conda environment specification
- `run.sh` - Setup script
- `run_pipeline.sh` - Automatic pipeline runner script
- `test_data/` - Test data (links to shared version)

## Notes

- Test data is referenced from the shared folder to maintain consistency across all platforms
- The nf-core RNA-seq pipeline uses a specific revision (3.14.0) to avoid version conflicts

---

For more details, see the shared README in `../../shared/nextflow/README.md`.
