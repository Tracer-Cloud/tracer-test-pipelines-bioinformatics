# Bioinformatics Pipeline - Bash Version

This directory contains a complete bioinformatics pipeline for variant calling using Docker containers.

## Pipeline Overview

The pipeline performs the following steps:

1. **FastQC** - Quality control on input FASTQ files
2. **TrimGalore** - Adapter trimming and quality filtering
3. **FastQC** - Quality control on trimmed files
4. **BWA-MEM** - Read alignment to reference genome
5. **SAMtools** - SAM to BAM conversion, sorting, and indexing
6. **GATK HaplotypeCaller** - Variant calling

## Prerequisites

- Docker installed and running
- wget (for downloading reference genome)
- At least 4GB of free disk space

## Quick Start

### 1. Setup Reference Genome

First, download and index the reference genome:

```bash
# Make scripts executable
chmod +x *.sh

# Setup reference genome (downloads ~60MB)
./setup_reference.sh
```

### 2. Generate Test Data (Optional)

If you don't have FASTQ files, generate small test files:

```bash
# Generate small test FASTQ files
./generate_test_data.sh
```

### 3. Run the Pipeline

#### Using existing test files:

```bash
# Run with existing test files
./run.sh small_test_fastq_1.fastq small_test_fastq_2.fastq output/test_output.vcf
```

#### Using generated test files:

```bash
# Run with generated test files
./run.sh test_data/test_1.fastq test_data/test_2.fastq output/test_output.vcf
```

#### Using your own FASTQ files:

```bash
# Run with your own FASTQ files
./run.sh /path/to/your/reads_1.fastq /path/to/your/reads_2.fastq /path/to/output/result.vcf
```

## Complete Setup and Run Commands

Here are all the commands you need to run the pipeline from scratch:

```bash
# 1. Navigate to the pipeline directory
cd pipelines/codespaces/bash

# 2. Make all scripts executable
chmod +x *.sh

# 3. Setup reference genome (downloads and indexes)
./setup_reference.sh

# 4. Generate test data (optional - skip if you have your own FASTQ files)
./generate_test_data.sh

# 5. Run the pipeline with existing test files
./run.sh small_test_fastq_1.fastq small_test_fastq_2.fastq output/test_output.vcf

# OR run with generated test files
./run.sh test_data/test_1.fastq test_data/test_2.fastq output/test_output.vcf
```

## Output Files

The pipeline generates several output files:

- **FastQC reports**: HTML and ZIP files with quality metrics
- **Trimmed FASTQ files**: `*_val_1.fq` and `*_val_2.fq`
- **Alignment files**: SAM, BAM, and indexed BAM files
- **VCF file**: Final variant calls

## Troubleshooting

### Common Issues:

1. **Docker not running**: Make sure Docker is installed and running
2. **Insufficient disk space**: The pipeline needs at least 4GB free space
3. **Reference files missing**: Run `./setup_reference.sh` first
4. **Permission denied**: Make scripts executable with `chmod +x *.sh`

### Check Docker Status:

```bash
docker --version
docker ps
```

### Check Available Space:

```bash
df -h
```

### Verify Reference Files:

```bash
ls -la reference/
```

## File Structure

```
pipelines/codespaces/bash/
├── run.sh                    # Main pipeline script
├── setup_reference.sh        # Reference genome setup
├── generate_test_data.sh     # Test data generator
├── small_test_fastq_1.fastq  # Existing test file 1
├── small_test_fastq_2.fastq  # Existing test file 2
├── reference/                # Reference genome files (created by setup)
│   ├── chr19.fa             # Reference genome
│   ├── chr19.fa.bwt         # BWA index
│   ├── chr19.fa.fai         # FASTA index
│   └── chr19.dict           # FASTA dictionary
├── test_data/               # Generated test files (optional)
│   ├── test_1.fastq
│   └── test_2.fastq
└── output/                  # Pipeline outputs (created automatically)
    ├── test_output.vcf      # Final variant calls
    ├── FastQC reports
    ├── BAM files
    └── Other intermediate files
```

## Performance Notes

- The pipeline uses chromosome 19 for testing (smaller, faster)
- For production use, consider using the full human genome
- Docker images are pulled automatically on first run
- Total runtime: ~10-30 minutes depending on data size and system performance
