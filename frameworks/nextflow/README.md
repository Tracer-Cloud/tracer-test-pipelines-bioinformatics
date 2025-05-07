


# Tracer Nextflow Pipelines

This directory contains reproducible pipeline wrappers for testing client workloads with [Nextflow](https://www.nextflow.io/) using the nf-core ecosystem.

## Folder Structure

- `tracer-test-pipelines-bioinformatics/frameworks/nextflow/pipelines/nf-core/*`: Cloned nf-core pipelines as Git submodules.
- `tracer-test-pipelines-bioinformatics/frameworks/nextflow/nextflow-config/`: Contains configuration files for local and AWS Batch execution.
- `tracer-test-pipelines-bioinformatics/Makefile`: Provides convenient targets to run pipeline tests locally or via AWS Batch.
- `tracer-test-pipelines-bioinformatics/spack.yaml`: Defines the environment and dependencies used to run pipelines with Spack.

## Setup

To initialize the environment and download dependencies using [Spack](https://spack.io/), run from the project root:

```bash
make setup_environment
```

This will:

* Download the specified Spack release.
* Initialize the environment.
* Install dependencies including Java and Nextflow.

## Running Pipelines

The Makefile contains the following targets for running test executions:

* **Sarek**

  * `make test_sarek`
  * `make test_sarek_aws_batch`
  * `make test_full_sarek_aws_batch`

* **RNA-seq**

  * `make test_rnaseq`
  * `make test_rnaseq_aws_batch`
  * `make test_full_rnaseq_aws_batch`

* **Proteinfold**

  * `make test_proteinfold`
  * `make test_proteinfold_aws_batch`
  * `make test_full_proteinfold_aws_batch`

Each target will:

* Activate the Spack environment.
* Run the corresponding Nextflow pipeline using configs from `nextflow-config/`.

## How to run on Tracer Sandbox

Only the test_rnaseq pipeline is tested to run completely on our sandbox

> Only the `test_rnaseq` pipeline has been verified to run successfully on the Tracer sandbox environment.

To run it:

```bash
make test_rnaseq
```

---

