# Minimal Nextflow Pipeline with Conda

A simple, reproducible Nextflow pipeline that runs on Mac using Conda instead of Docker containers.

## What This Pipeline Does

1. **FastQC Analysis** - Runs quality control on FASTA files
2. **Sequence Counting** - Counts sequences in each FASTA file  
3. **Summary Report** - Generates a summary of all processed files

## Prerequisites

- **Conda/Miniconda** - Install via Homebrew: `brew install miniconda`
- **Nextflow** - Will be installed automatically by setup script

## Quick Start

1. **Clone or download these files to a directory**

2. **Run the setup script:**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```

3. **Run the pipeline:**
   ```bash
   nextflow run main.nf --input 'test_data/*.fasta'
   ```

## Pipeline Structure

```
.
├── main.nf              # Main pipeline workflow
├── nextflow.config      # Configuration and parameters
├── environment.yml      # Conda environment specification
├── setup.sh            # Automated setup script
├── test_data/          # Sample FASTA files (created by setup)
└── results/            # Output directory
    ├── fastqc/         # FastQC HTML reports
    ├── counts/         # Sequence count files
    ├── summary_report.txt
    ├── pipeline_report.html
    └── timeline.html
```

## Usage Examples

### Basic Usage
```bash
nextflow run main.nf --input 'test_data/*.fasta'
```

### Mac Optimized (more CPU/memory)
```bash
nextflow run main.nf --input 'test_data/*.fasta' -profile mac
```

### Custom Output Directory
```bash
nextflow run main.nf --input 'test_data/*.fasta' --outdir my_results
```

### Debug Mode (keeps work files)
```bash
nextflow run main.nf --input 'test_data/*.fasta' -profile debug
```

### Resume Failed Pipeline
```bash
nextflow run main.nf --input 'test_data/*.fasta' -resume
```

## Configuration Profiles

- **standard** (default) - 2 CPUs, 4GB RAM
- **mac** - 4 CPUs, 8GB RAM, longer timeouts
- **debug** - Minimal resources, keeps work files

## Parameters

- `--input` - Input FASTA files pattern (default: `*.fasta`)
- `--outdir` - Output directory (default: `results`)

## Output Files

- **fastqc/** - HTML quality control reports
- **counts/** - Text files with sequence counts
- **summary_report.txt** - Combined summary of all files
- **pipeline_report.html** - Nextflow execution report
- **timeline.html** - Pipeline execution timeline

## Conda Environment

The pipeline automatically creates a Conda environment with:
- FastQC 0.12.1
- Core Unix utilities
- SeqKit (for advanced FASTA handling)
- Python 3.9 + BioPython

## Troubleshooting

### Environment Creation Takes Long
- The first run creates the Conda environment (can take 10-30 minutes)
- Subsequent runs reuse the environment and are much faster
- Using Mamba (installed by setup script) speeds this up

### Pipeline Fails
- Check the `.nextflow.log` file for detailed errors
- Run with `-profile debug` to keep work files for inspection
- Use `-resume` to restart from the last successful step

### No Input Files Found
- Check your file pattern: `--input 'path/to/*.fasta'`
- Ensure FASTA files exist in the specified location
- Use absolute paths if relative paths don't work

## Extending the Pipeline

To add new processes:

1. Add process definition to `main.nf`
2. Update the workflow section to include the new process
3. Add any new software requirements to `environment.yml`
4. Test with `-profile debug` first

## Benefits of This Approach

✅ **No Docker complexity** - Pure Conda-based approach  
✅ **Reproducible** - Pinned software versions  
✅ **Mac-friendly** - Optimized for macOS development  
✅ **Portable** - Works on any system with Conda  
✅ **Fast** - Local execution with Conda environments  
✅ **Debuggable** - Easy to inspect and modify