
# nf-core/sarek

This directory runs the [`nf-core/sarek`](https://github.com/nf-core/sarek) pipeline using a lightweight local setup powered by [Spack](https://spack.io/) and Nextflow.

## Overview

The `sarek` pipeline is designed for germline or somatic variant calling from raw sequencing data. This wrapper provides a minimal test profile to validate the integration inside Tracer environments.

## Requirements

- Docker (required for pipeline execution)
- Nextflow (installed via Spack)
- Spack 0.23.0 (downloaded automatically)

## Setup

To bootstrap the environment and install dependencies:

```bash
make setup_environment
````

This will:

* Download and extract Spack
* Bootstrap Spack's compiler/runtime support
* Install required dependencies in-place

## Running the Pipeline

To execute the Sarek pipeline using a lightweight test profile:

```bash
make test_sarek
```

This will:

* Activate the Spack environment
* Run `nf-core/sarek` via Nextflow
* Use parameters from `../config/sarek-params.json`
* Use `../config/local.config` as configuration
* Enable the `docker,arm,test` profile combination

