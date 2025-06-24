# Nextflow Pipeline with Pixi Environment Management

This directory contains a minimal Nextflow pipeline that uses [Pixi](https://pixi.sh) for dependency management. Pixi provides faster, more reliable environment management with better reproducibility compared to Conda.

## Quick Start

Make sure the tracer daemon is running before executing the pipeline so that the tools are recognized.

### Prerequisites

- [Pixi](https://pixi.sh) installed on your system
- macOS ARM64 (Apple Silicon)

### Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Run the Pipeline

```bash
pixi install

pixi run pipeline

pixi run pipeline-custom "path/to/input/*.fasta" "path/to/output"
```

## 🧬 Pipeline Description

This pipeline processes FASTA files using FastQC and STAR, performing the following steps:

1. **Version Check**: Displays versions of installed tools (FastQC and STAR)
2. **FASTA Statistics**:
   - Counts sequences in each file
   - Calculates total sequence length
   - Computes average sequence length
   - Lists sequence IDs
3. **Sequence Counting**: Generates separate count files for each FASTA input
