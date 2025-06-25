# Multi-step Nextflow Pipeline with Conda (macOS ARM64)

## Quick Start

**Make sure the tracer daemon is running before executing the pipeline so that the tools are recognized.**

### Check Tracer Status

- To monitor the Tracer daemon and running processes, open the Tracer dashboard in your browser (usually at http://localhost:3000).
- Alternatively, you can view active processes and tracer status in your terminal:

```bash
tracer info
```

### Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) installed on your system
- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- macOS ARM64 (Apple Silicon)

### Install Dependencies

```bash
# Create and activate the Conda environment
conda env create -f environment.yml
conda activate nextflow-minimal
```

### Run the Pipeline

Smaller pipeline

```bash
nextflow run main.nf
```

Rnaseq pipeline

```bash
mkdir -p rnaseq_results && cd rnaseq_results && nextflow run https://github.com/nf-core/rnaseq -c ../custom.config -profile test --outdir results -resume
```

> After running the pipeline, you can view the tools and processes in the Tracer dashboard or by running `tracer info` in your terminal.
