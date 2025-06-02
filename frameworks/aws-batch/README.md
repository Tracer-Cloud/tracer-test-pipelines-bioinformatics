# Tracer AWS Batch Pipelines

This directory enables distributed, cloud-based testing of nf-core bioinformatics pipelines using [Nextflow](https://www.nextflow.io/) with [AWS Batch](https://aws.amazon.com/batch/).

## Folder Structure

* `frameworks/nextflow/aws-batch/`: Contains this Makefile and batch-specific test targets.
* `frameworks/config/nextflow/`: Contains AWS Batch configuration (`batch.config`, GPU config, params).
* `frameworks/pipelines/nf-core/*`: Cloned nf-core pipelines (e.g. `rnaseq`, `sarek`, `proteinfold`).
* `frameworks/spack.yaml`: Declares Spack dependencies (Java, Nextflow, etc).

## Setup

Before running any pipeline, set up the environment from this directory:

```bash
make setup_environment
```

This will:

* Download Spack and bootstrap it.
* Install required tools.
* Prepare the Nextflow environment.

## Running Pipelines on AWS Batch

Each test target submits a specific pipeline run to AWS Batch using a profile in `batch.config`. The pipelines currently supported:

| Pipeline           | Target                                 | Description                         |
| ------------------ | -------------------------------------- | ----------------------------------- |
| Sarek              | `make test_sarek_aws_batch`            | Run Sarek pipeline (minimal test)   |
| Sarek (full)       | `make test_full_sarek_aws_batch`       | Run extended Sarek test             |
| RNA-seq            | `make test_rnaseq_aws_batch`           | Run RNA-seq pipeline (minimal test) |
| RNA-seq (full)     | `make test_full_rnaseq_aws_batch`      | Run extended RNA-seq test           |
| Proteinfold        | `make test_proteinfold_aws_batch`      | Run GPU-enabled Proteinfold test    |
| Proteinfold (full) | `make test_full_proteinfold_aws_batch` | Full GPU Proteinfold run            |

## How It Works

* **Nextflow** orchestrates the pipeline and submits jobs to AWS Batch.
* **CloudFormation** provisions compute environments with appropriate IAM roles, queues, and resources.
* **Tracer Agent** is auto-installed on each EC2 instance via bootstrap scripts.
* **Session Tracking** is enabled via `workflow.sessionId` passed as an environment variable.
* **GPU Pipelines** (e.g. Proteinfold) include `proteinfold.config` for resource definitions.

## Execution Environments

### Option 1: Your AWS Account

1. Clone this repo and enter the `aws-batch` directory.
2. Make sure AWS credentials are configured.
3. Run a test pipeline like:

```bash
make test_rnaseq_aws_batch
```

View pipeline metrics and logs in your Tracer observability dashboard.

### Option 2: Tracer Sandbox (Preconfigured)

1. Spin up an EC2 Sandbox environment using the provided EC2 launch template.
2. Access the sandbox and switch to the correct user directory:

```bash

sudo su - ubuntu && cd tracer-test-pipelines-bioinformatics

# Initialize the pipeline
tracer init --pipeline-name aws_batch_test \
    --environment sandbox \
    --user-operator vincent \
    --pipeline-type aws_batch_rnaseq

# Run the test
make test_rnaseq_aws_batch
```

---
