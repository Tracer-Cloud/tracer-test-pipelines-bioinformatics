


# Tracer Nextflow Pipelines

This directory contains reproducible pipeline wrappers for testing client workloads with [Nextflow](https://www.nextflow.io/) using the nf-core ecosystem.

## Folder Structure

* `integrations/pipelines/nf-core/*`: Cloned nf-core pipelines (e.g., `sarek`, `rnaseq`, `proteinfold`) as Git submodules.
* `integrations/config/nextflow/`: Contains configuration files and parameter sets for local execution.
* `integrations/nextflow/`: Contains per-pipeline Makefiles for development and testing.
* `integrations/spack.yaml`: Defines Spack environment dependencies, including Java and Nextflow.
* `integrations/nextflow/Makefile`: Contains targets to automate testing and execution of pipelines.

## Setup

To initialize the environment and install required tools via [Spack](https://spack.io/), run from within the `integrations/nextflow/` directory:

```bash
make setup_environment
```

This will:

* Download the specified Spack release (version `0.23.0`).
* Extract and bootstrap Spack.
* Install dependencies such as Java and Nextflow.
* Activate the environment in-place.

## Running Pipelines

Each pipeline has a Makefile target that:

* Activates the local Spack environment.
* Runs the appropriate nf-core pipeline using `nextflow`.
* Applies parameters and configuration from `integrations/config/nextflow/`.

Available targets:

| Pipeline       | Target Command          | Description                      |
| -------------- | ----------------------- | -------------------------------- |
| Sarek          | `make test_sarek`       | Run minimal test profile locally |
| RNA-seq        | `make test_rnaseq`      | Run minimal test profile locally |
| RNA-seq (full) | `make test_rnaseq_full` | Run extended test locally        |
| Proteinfold    | `make test_proteinfold` | Run minimal test profile locally |

Each target references:

* `local.config` from `config/nextflow`
* A `*-params.json` file containing pipeline-specific inputs

## Example: Running on Tracer Sandbox

To test the RNA-seq pipeline on the Tracer sandbox:

```bash
cd integrations/nextflow
make test_rnaseq
```

*Note: Only `test_rnaseq` has been verified end-to-end in the sandbox environment so far.*

---
