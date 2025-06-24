# Linux x86-64 Ubuntu Nextflow Pipeline with Docker

This directory contains a Nextflow pipeline setup for Linux x86-64 Ubuntu systems using Docker for containerized execution.

## Quick Start

### Prerequisites
- Docker installed and running
- Docker daemon accessible to your user

### Build the Docker Image

```bash
# Build the custom image with all bioinformatics tools
docker build -t rnaseq-tools:latest .
```

### Run nf-core Pipelines

```bash
# Run the nf-core rnaseq pipeline with Docker
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results
```

## Files

- `dockerfile` - Docker image with all required bioinformatics tools
- `custom.config` - Nextflow configuration optimized for Docker execution
- `README.md` - This documentation

## Docker Image Contents

The Docker image includes:
- **Ubuntu 22.04** base system
- **Miniconda** for package management
- **Bioinformatics tools**: FastQC, STAR, HISAT2, Salmon, Samtools, BBMap, Trim-galore, etc.
- **Nextflow** workflow manager
- **System utilities**: gawk, bc, curl, wget

## Troubleshooting

### Docker Permission Issues
The custom.config includes Docker settings to handle file ownership:
- `fixOwnership = true`
- `runOptions = '--user $(id -u):$(id -g)'`

### Resource Constraints
The configuration limits resources to work on smaller systems:
- Max 2 CPUs per process
- Max 3 GB memory per process
- Specific limits for memory-intensive processes

### Docker Not Running
Make sure Docker is installed and running:
```bash
sudo systemctl start docker
sudo usermod -aG docker $USER
```

## Usage Example

```bash
# Build the image
docker build -t rnaseq-tools:latest .

# Run a test pipeline
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results

# Check results
ls -la results/
```
