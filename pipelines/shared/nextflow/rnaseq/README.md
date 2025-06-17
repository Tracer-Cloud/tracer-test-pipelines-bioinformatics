# nf-core RNA-seq Pipeline

This directory contains a reproducible setup for running the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline using [Nextflow](https://www.nextflow.io/) with Tracer integration.

## Overview

The RNA-seq pipeline processes raw RNA sequencing data and includes quality control, alignment, and quantification steps. This setup uses a shared Spack environment and central configuration files for consistency across environments.

## Requirements

- Spack environment (run `make setup_environment` from the `nextflow/` root)
- Configuration and parameter files are expected at `../config/`
- Docker or Singularity (depending on the selected profile)

## Usage

Run from within this directory:

```bash
make test_rnaseq
````

This executes the pipeline with:

* `../config/local.config` as the config
* `../config/rnaseq-params.json` as the input parameter file
* `-profile test`

To run the extended/full test:

```bash
make test_rnaseq_full
```

This switches the profile to `test_full`.

## Makefile Targets

| Target                  | Description                                |
| ----------------------- | ------------------------------------------ |
| `make test_rnaseq`      | Runs a minimal test pipeline               |
| `make test_rnaseq_full` | Runs extended version of the test pipeline |





