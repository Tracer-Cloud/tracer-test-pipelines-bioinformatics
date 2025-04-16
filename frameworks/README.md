# Tracer Workload Validation & Installation Guide
This document provides an overview of Tracer workloads, tracks technical development progress, and ensures validated workflow coverage.

## Workloads Validation Status
| Workload                | Status           | Example Configuration |
|-------------------------|------------------|------------------------|
| AWS Batch               | Validated        | [Link](./airflow/README.md) |
| Bash (RNA-seq, ChIP-seq)| Validated        | [Link](./bash/README.md)    |
| Nextflow on EC2         | Validated        | [Link](./nextflow)          |
| Airflow                 | Validated        | [Link](./airflow/README.md) |
| Slurm                   | Not Tested       | [Link](./slurm/README.md)   |
| R Bioconductor          | Not Tested       | [Link](#)                   |
| AlphaFold               | Not Tested       | [Link](#)                   |
| OpenFold                | Not Tested       | [Link](#)                   |

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
