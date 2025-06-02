# nf-core/proteinfold

This directory wraps the [`nf-core/proteinfold`](https://github.com/nf-core/proteinfold) pipeline using Spack and Nextflow for simple reproducibility in local and sandbox environments.

## Overview

The `proteinfold` pipeline uses deep learning to predict the 3D structure of proteins from amino acid sequences. This wrapper ensures a minimal test profile is runnable in containerized setups like Tracer.

## Requirements

- Docker (required for running the containerized pipeline)
- Nextflow (installed via Spack)
- Spack 0.23.0 (bootstrapped automatically)

## Setup

Run the following to initialize your local Spack environment:

```bash
make setup_environment
````

This:

* Downloads and configures Spack
* Activates the project environment
* Installs Java and Nextflow (and related tools)

## Running the Pipeline

To execute a minimal test profile:

```bash
make test_proteinfold
```

This will:

* Activate the Spack environment
* Use `../config/proteinfold-params.json` as input
* Apply `../config/local.config` as config
* Run `nf-core/proteinfold` with the `docker,arm,test` profile


