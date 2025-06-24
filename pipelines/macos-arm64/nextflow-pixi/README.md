# Nextflow Pipeline with Pixi Environment Management

This directory contains a minimal Nextflow pipeline that uses [Pixi](https://pixi.sh) for dependency management. Pixi provides faster, more reliable environment management with better reproducibility compared to Conda.

## Quick Start

Make sure the tracer daemon is running before executing the pipeline so that the tools are recognized.

### Prerequisites

- macOS ARM64 (Apple Silicon)
- [Pixi](https://pixi.sh) installed on your system
- Docker

### Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Run the Pipeline

```bash
pixi install

pixi run pipeline

pixi run rnaseq

pixi run pipeline-custom "path/to/input/*.fasta" "path/to/output"
```

> After running the pipeline, you can view the tools and processes in the Tracer dashboard or by running `tracer info` in your terminal.
