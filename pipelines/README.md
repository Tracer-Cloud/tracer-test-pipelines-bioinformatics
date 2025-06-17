# Tracer Workload Validation & Installation Guide

This document provides an overview of Tracer workloads, tracks technical development progress, and ensures validated workflow coverage.

## Workloads Validation Status

| Workload                 | Status     | Example Configuration                |
| ------------------------ | ---------- | ------------------------------------ |
| AWS Batch                | Validated  | [Link](./aws-batch/README.md)        |
| Bash (RNA-seq, ChIP-seq) | Validated  | [Link](./bash/README.md)             |
| Nextflow                 | Validated  | [Link](./nextflow)                   |
| Nextflow (Conda)         | Validated  | [Link](./nextflow-conda-macos-arm64) |
| Nextflow (AWS x86)       | Validated  | [Link](./nextlow-aws-x86-ubuntu)     |
| Airflow                  | Validated  | [Link](./airflow/README.md)          |
| CWL                      | Validated  | [Link](./cwl)                        |
| WDL                      | Validated  | [Link](./wdl)                        |
| Slurm                    | Not Tested | [Link](./slurm/README.md)            |
| R Bioconductor           | Not Tested | [Link](#)                            |
| AlphaFold                | Not Tested | [Link](#)                            |
| OpenFold                 | Not Tested | [Link](#)                            |

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
