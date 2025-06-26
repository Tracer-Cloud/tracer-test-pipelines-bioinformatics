# Codespaces Pipeline Testing

This directory contains Nextflow pipeline configurations optimized for GitHub Codespaces environments, supporting both Pixi and Conda dependency management systems.

## Overview

The codespaces setup provides two different approaches to running Nextflow pipelines:

1. **Pixi Environment** (`nextflow-pixi/`) - Uses Pixi for dependency management
2. **Conda Environment** (`nextflow-conda/`) - Uses Conda for dependency management

Both environments run the same version check workflow to verify bioinformatics tool availability.

## Quick Start

### Option 1: Pixi Environment (Recommended)

```bash
cd pipelines/codespaces/nextflow-pixi
bash run.sh
```

### Option 2: Conda Environment

```bash
cd pipelines/codespaces/nextflow-conda
bash run.sh
```

## What Each Environment Does

### Pixi Environment (`nextflow-pixi/`)

- **Dependency Management**: Uses Pixi for fast, reproducible dependency resolution
- **Tools Included**: Nextflow, Java, FastQC, STAR, Samtools, BWA
- **Memory Optimized**: Configured for 800MB memory usage
- **Features**:
  - Automatic Pixi installation if not present
  - Environment verification with `pixi run check-env`
  - Pipeline execution with `pixi run pipeline`
  - Cleanup with `pixi run clean`

### Conda Environment (`nextflow-conda/`)

- **Dependency Management**: Uses Conda/Mamba for traditional package management
- **Tools Included**: Nextflow, Java, FastQC, STAR, Samtools
- **Memory Optimized**: Configured for 800MB memory usage
- **Features**:
  - Smart conda detection (system conda, miniconda, etc.)
  - Automatic environment creation/update
  - Tool verification before pipeline execution
  - Robust error handling

## Expected Output

Both environments will generate:

```
results/
├── tool_versions.txt          # Version information for all tools
└── pipeline_report.html       # Nextflow execution report

logs/
├── nextflow.log              # Detailed execution log
├── timeline_report.html      # Timeline of process execution
└── trace.txt                 # Execution trace
```

## GitHub Actions Integration

The `.github/workflows/codespaces.yml` workflow automatically tests both environments:

- **Triggers**: Push to `pipelines/codespaces/**`, manual dispatch, daily schedule
- **Tests**: Both pixi and conda environments
- **Artifacts**: Results and logs are uploaded as artifacts
- **Cleanup**: Automatic cleanup after testing

## Troubleshooting

### Memory Issues

If you encounter memory errors, the configurations are already optimized for low-memory environments (800MB). If issues persist:

1. Check available memory: `free -h`
2. Ensure no other processes are consuming memory
3. Consider using a larger codespace instance

### Tool Not Found Errors

If tools are not found:

1. **Pixi**: Run `pixi run check-env` to verify all tools
2. **Conda**: Check if conda environment is activated: `conda info --envs`

### Pipeline Failures

If the pipeline fails:

1. Check logs in the `logs/` directory
2. Verify input data exists in `test_data/`
3. Ensure all dependencies are properly installed

## Environment Comparison

| Feature               | Pixi            | Conda           |
| --------------------- | --------------- | --------------- |
| Installation Speed    | ⚡ Fast         | 🐌 Slower       |
| Dependency Resolution | 🔄 Reproducible | 🔄 Reproducible |
| Memory Usage          | 💾 Low          | 💾 Low          |
| Tool Availability     | ✅ All tools    | ✅ All tools    |
| Error Handling        | 🛡️ Robust       | 🛡️ Robust       |
| Cleanup               | 🧹 Automatic    | 🧹 Manual       |

## Development

To modify the pipelines:

1. **Pixi**: Edit `pixi.toml` for dependencies, `main.nf` for workflow
2. **Conda**: Edit `environment.yml` for dependencies, `main.nf` for workflow
3. **Configuration**: Edit `nextflow.config` for execution settings

## Support

For issues or questions:

- Check the main [README.md](../../README.md)
- Review the workflow logs in GitHub Actions
- Ensure your codespace has sufficient resources
