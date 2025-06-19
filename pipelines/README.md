# Tracer Bioinformatics Test Pipelines

This repository contains validated test pipelines for various bioinformatics platforms and environments, with **Pixi as the default dependency management solution**.

## 🚀 Quick Start (Recommended)

### Pixi-based Pipelines (Default)
```bash
# For Apple Silicon Macs
cd macos-arm64/nextflow-pixi && pixi run pipeline

# For Intel Macs
cd macos-intel-x86/nextflow-pixi && pixi run pipeline
```

## Workloads Validation Status

| Workload                 | Status     | Dependency Manager | Example Configuration                |
| ------------------------ | ---------- | ------------------ | ------------------------------------ |
| **Nextflow (Pixi ARM64)** | ✅ **Default** | Pixi | [Link](./macos-arm64/nextflow-pixi) |
| **Nextflow (Pixi Intel)** | ✅ **Default** | Pixi | [Link](./macos-intel-x86/nextflow-pixi) |
| AWS Batch                | ✅ Validated  | Direct | [Link](./aws-batch/README.md)        |
| Bash (RNA-seq, ChIP-seq) | ✅ Validated  | System | [Link](./bash/README.md)             |
| Nextflow (Config)        | ✅ Validated  | Various | [Link](./nextflow)                   |
| Airflow                  | ✅ Validated  | Conda | [Link](./airflow/README.md)          |
| CWL                      | ✅ Validated  | System | [Link](./shared/cwl)                        |
| WDL                      | ✅ Validated  | System | [Link](./wdl)                        |
| Slurm                    | ⏳ Not Tested | System | [Link](./slurm/README.md)            |
| R Bioconductor           | ⏳ Not Tested | R | [Link](#)                            |
| AlphaFold                | ⏳ Not Tested | Conda/Docker | [Link](#)                            |
| OpenFold                 | ⏳ Not Tested | Conda/Docker | [Link](#)                            |

## ⚡ Why Pixi? (Migration from Conda)

We've migrated from Conda to Pixi as our default dependency manager:

### **Performance Benefits:**
- **5-10x faster** environment creation (30-60s vs 2-5 minutes)
- **Better dependency resolution** with fewer conflicts
- **Faster CI/CD** with reliable caching

### **Developer Experience:**
- **Task-based workflow**: `pixi run pipeline`, `pixi run test`
- **Automatic environment activation**: No manual conda activate
- **Built-in lock files**: Guaranteed reproducibility

### **Migration Path:**
- ❌ **Old**: `pipelines/macos-*/nextflow-conda` (removed)
- ✅ **New**: `pipelines/macos-*/nextflow-pixi` (default)

# Core Functionality Requirements

- Identification of individual pipeline runs and process stages
- Core metrics tracking:
  - Execution duration
  - CPU, Memory, and Disk usage
  - Cost estimation

## Link To Repository With Examples:

- https://github.com/Tracer-Cloud/tracer-workflow-templates

# Instructions

## Instruction: AWS Batch With Nextflow

- Needs specification

## Instruction: Bash Script

- https://github.com/Tracer-Cloud/tracer-workflow-templates

## Instruction: Airflow

- Needs specification

## Instruction: Slurm

- Needs specification.
