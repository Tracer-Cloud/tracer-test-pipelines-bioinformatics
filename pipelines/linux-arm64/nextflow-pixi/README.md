# Nextflow Pipeline with Pixi Environment Management (Linux ARM64)

This directory contains a minimal Nextflow pipeline that uses [Pixi](https://pixi.sh) for dependency management on Linux ARM64 systems. Pixi provides faster, more reliable environment management with better reproducibility.

## 🚀 Quick Start

### Prerequisites
- [Pixi](https://pixi.sh) installed on your system
- Linux ARM64 (aarch64) architecture

### Install Pixi
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Run the Pipeline
```bash
# Navigate to this directory
cd pipelines/linux-arm64/nextflow-pixi

# Install dependencies and run pipeline
pixi run pipeline

# Or run with custom parameters
pixi run pipeline-custom --input "custom_data/*.fasta" --outdir "custom_results"
```

## 📁 Directory Structure

```
nextflow-pixi/
├── pixi.toml              # Pixi project configuration
├── nextflow.config        # Nextflow configuration
├── main.nf               # Main Nextflow script
├── run.sh                # Pipeline runner script
├── test_data/            # Sample FASTA files (created by setup)
│   ├── sample1.fasta
│   └── sample2.fasta
└── README.md             # This file
```

## 🔧 Available Tasks

Pixi provides several predefined tasks:

```bash
# Run the main pipeline
pixi run pipeline

# Run with custom parameters
pixi run pipeline-custom

# Run tests
pixi run test

# Clean up results
pixi run clean

# Check environment
pixi run check-env

# Development workflow (clean + setup + check + test)
pixi run dev
```

## 🧬 Pipeline Description

This pipeline processes FASTA files and performs:

1. **Version Check**: Displays versions of FastQC and STAR
2. **FastQC Analysis**: Runs FastQC on input files (with FASTA adaptation)
3. **STAR Information**: Shows STAR version and capabilities
4. **FASTA Statistics**: Calculates sequence counts, lengths, and statistics

### Outputs
- `results/versions/`: Tool version information
- `results/fastqc/`: FastQC analysis results and summaries
- `results/star/`: STAR version and information
- `results/fasta_stats/`: Detailed statistics for each FASTA file
- `logs/`: Pipeline execution logs and reports

## ⚡ Advantages of Pixi over Conda

| Feature | Conda | Pixi |
|---------|-------|------|
| Speed | Slow dependency resolution | Fast dependency resolution |
| Lock files | Manual (conda-lock) | Automatic |
| Cross-platform | Good | Excellent |
| Reproducibility | Good | Excellent |
| ARM64 support | Limited | Native |

## 🏗️ Architecture-Specific Features

This Linux ARM64 version includes:

- **Native ARM64 binaries**: All tools compiled for aarch64
- **Optimized resource allocation**: ARM64-specific CPU and memory settings
- **Platform detection**: Automatic architecture detection and reporting
- **ARM64 profile**: Specialized execution profile for ARM64 systems

## 🔧 Configuration

The pipeline supports multiple execution profiles:

- `standard`: Default local execution
- `debug`: Verbose output with detailed logging
- `fast`: Higher resource allocation
- `arm64`: ARM64-optimized settings

## 🧪 Testing

```bash
# Run basic tests
pixi run test

# Run with debug output
pixi run pipeline --profile debug

# Check environment
pixi run check-env
```

## 🔗 Related Files

- `../linux-intel-x86/nextflow-pixi/`: Intel x86 Linux version
- `../../macos-arm64/nextflow-pixi/`: macOS ARM64 version
- `../../codespaces/bash/`: Codespaces Pixi runner

## 📊 Performance Notes

On Linux ARM64 systems:
- FastQC: ~30 seconds for typical FASTA files
- STAR: Version check and info generation ~5 seconds
- Overall pipeline: ~1-2 minutes for test data

## 🐛 Troubleshooting

Common issues and solutions:

1. **Pixi not found**: Ensure Pixi is installed and in PATH
2. **Architecture mismatch**: Verify you're on Linux ARM64 (aarch64)
3. **Permission errors**: Ensure write permissions in output directory
4. **Memory issues**: Use `--profile fast` for more resources

## 🚀 Next Steps

- Extend with additional bioinformatics tools
- Add genome indexing for STAR
- Implement quality control workflows
- Add support for paired-end data
