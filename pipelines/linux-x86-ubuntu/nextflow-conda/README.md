# Linux x86-64 Ubuntu Nextflow Pipeline with Conda/Pixi

This directory contains a Nextflow pipeline setup for Linux x86-64 Ubuntu systems supporting both Conda and Pixi package management for fast, reliable dependency management.

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
pixi run pipeline-custom my_results

# Clean up
pixi run clean
```

## Alternative: Automated Setup (Bash Script)

The `run.sh` script now installs both Pixi and Conda automatically:

1. **Run the setup script (installs Pixi, Conda, Java, and Nextflow):**
   ```bash
   ./run.sh
   ```

2. **After setup, you can use either approach:**

   **Option A: Use Pixi (recommended):**
   ```bash
   pixi install
   pixi run pipeline
   ```

   **Option B: Use Conda:**
   ```bash
   source ~/.bashrc
   conda --version
   nextflow run main.nf --outdir results
   ```

   **Option C: Use the conda activation helper:**
   ```bash
   source activate-conda.sh
   ```

## Files

- `pixi.toml` - Pixi configuration with dependencies and tasks (recommended)
- `main.nf` - Nextflow pipeline definition (simple tool version checks)
- `nextflow.config` - Nextflow configuration
- `custom.config` - Custom configuration for nf-core pipelines
- `test_data/` - Sample test data for the pipeline
- `run.sh` - Automated setup script that installs Pixi, Conda, Java, and Nextflow
- `activate-conda.sh` - Helper script to activate conda in current shell session

## Troubleshooting

### Conda command not found after running run.sh

This is expected behavior. The `run.sh` script installs conda and configures it, but you need to activate it in your shell session:

**Option 1:** Source your bashrc file
```bash
source ~/.bashrc
```

**Option 2:** Use the activation helper script
```bash
source activate-conda.sh
```

**Option 3:** Restart your terminal/shell session

### Why does this happen?

The setup script modifies your `~/.bashrc` file to include conda in your PATH, but these changes only take effect in new shell sessions or when you explicitly source the bashrc file.

## Usage

### Pixi Tasks (Recommended)

```bash
# Check all available tasks
pixi task list

# Setup environment and directories
pixi run setup

# Check tool versions and environment
pixi run check-env

# Run the main pipeline
pixi run pipeline

# Run pipeline with custom output
pixi run pipeline-custom my_custom_results

# Run development workflow (clean + setup + check + test)
pixi run dev

# Run nf-core RNA-seq pipeline
pixi run nf-core-rnaseq

# Clean up all generated files
pixi run clean
```

### Pipeline Description

This directory contains a simple Nextflow pipeline that checks versions of common bioinformatics tools:

The pipeline will:
- Check versions of FastQC, STAR, Samtools, BWA, and GATK
- Create a `tool_versions.txt` file in the results directory
- Continue even if some tools are not available (uses `errorStrategy = 'ignore'`)
- Use 1GB memory per process step

### Manual Nextflow Execution

If you prefer to run nextflow directly (after `pixi install`):

```bash
# Run the local pipeline
nextflow run main.nf --outdir results

# Run nf-core rnaseq pipeline
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results -resume
```

### Legacy Conda Usage

If using the bash setup instead of pixi:

```bash
conda --version
conda list
conda create -n myenv python=3.9
conda activate myenv
```

The `custom.config` file is included in this directory and contains:
- `process.errorStrategy = 'ignore'` - Allows the pipeline to continue even if some processes fail

## Why Pixi?

Pixi offers several advantages over traditional conda/bash setup:

- **Faster**: Parallel dependency resolution and installation
- **Reproducible**: Lock files ensure exact dependency versions
- **Isolated**: Each project has its own environment
- **Cross-platform**: Works consistently across Linux, macOS, and Windows
- **Task Management**: Built-in task runner for common workflows
- **No Activation**: Tools are automatically available when running tasks

## Dependencies

The pixi environment includes:
- **Nextflow** (>=25.4.4): Workflow management system
- **FastQC** (0.12.1): Quality control for sequencing data
- **STAR** (2.7.11b): RNA-seq aligner
- **Samtools** (1.21): SAM/BAM file manipulation
- **BWA** (0.7.18): DNA sequence aligner
- **GATK4** (4.6.1.0): Genome analysis toolkit
- **OpenJDK** (17.0.11): Java runtime for Nextflow and GATK