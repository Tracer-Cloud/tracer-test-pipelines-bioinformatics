# Nextflow Pipeline with Pixi Environment Management

This directory contains a minimal Nextflow pipeline that uses [Pixi](https://pixi.sh) for dependency management instead of Conda. Pixi provides faster, more reliable environment management with better reproducibility.

## 🚀 Quick Start

### Prerequisites
- [Pixi](https://pixi.sh) installed on your system
- macOS ARM64 (Apple Silicon)

### Install Pixi
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Run the Pipeline
```bash
# Navigate to this directory
cd pipelines/macos-arm64/nextflow-pixi

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
├── test_data/            # Sample FASTA files
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
2. **FASTA Statistics**: Calculates sequence counts, lengths, and statistics
3. **Sequence Counting**: Counts sequences in each FASTA file

### Outputs
- `results/fasta_stats/`: Detailed statistics for each FASTA file
- `results/counts/`: Sequence counts for each file
- `logs/`: Pipeline execution logs and reports

## ⚡ Advantages of Pixi over Conda

1. **Faster**: Rust-based solver is significantly faster than Conda
2. **More Reliable**: Better dependency resolution and conflict handling
3. **Reproducible**: Lock files ensure exact dependency versions
4. **Project-Centric**: Automatic environment activation per project
5. **Better CI/CD**: More predictable behavior in automated environments

## 🔄 Migration from Conda

This pipeline is equivalent to the `nextflow-conda` version but uses Pixi for dependency management:

### Key Differences:
- **No `environment.yml`**: Dependencies defined in `pixi.toml`
- **No conda configuration**: Nextflow config simplified
- **Task-based workflow**: Pixi tasks replace shell scripts
- **Automatic environment**: No manual activation needed

### Dependencies:
- FastQC 0.12.1
- STAR 2.7.11b
- coreutils 9.1
- bc (calculator)

## 🛠️ Development

### Adding Dependencies
```bash
# Add a new bioinformatics tool
pixi add samtools

# Add a development dependency
pixi add --feature dev jupyter
```

### Custom Tasks
Edit `pixi.toml` to add new tasks:
```toml
[tasks]
my-task = "echo 'Custom task'"
```

## 🐛 Troubleshooting

### Environment Issues
```bash
# Check environment status
pixi info

# Reinstall environment
pixi install --force

# Check tool availability
pixi run check-env
```

### Pipeline Issues
```bash
# Clean and restart
pixi run clean
pixi run dev

# Check logs
cat logs/nextflow.log
```

## 📊 Performance Comparison

| Aspect | Conda | Pixi |
|--------|-------|------|
| Environment creation | ~2-5 min | ~30-60 sec |
| Dependency resolution | Sometimes fails | More reliable |
| Lock files | conda-lock (separate tool) | Built-in |
| Cross-platform | Good | Excellent |
| CI/CD performance | Slow | Fast |

## 🔗 Related Files

- `../nextflow-conda/`: Original Conda-based version
- `../../workflows/macos-intel-x86.yml`: GitHub Actions workflow
- `../../workflows/mac-arm64.yml`: GitHub Actions workflow (ARM64)
