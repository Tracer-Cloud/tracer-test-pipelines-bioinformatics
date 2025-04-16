# Development Troubleshooting For AWS Batch 

This guide outlines methods for local and cloud-based testing of Tracer pipelines.

## Testing Locally (Cloud Services)
To run basic tests locally using cloud services, execute:

```bash
make test-tracer
```

## Testing with AWS Batch (This needs to be moved to the testing repository)
- AWS Batch is used to run scalable, distributed bioinformatics pipelines powered by Nextflow. 
- Here's a breakdown of how Tracer integrates with this setup for observability and performance monitoring:

**Nextflow Orchestrates the Workload:**
Nextflow is responsible for generating AWS Batch jobs that distribute different pipeline processes across multiple EC2 instances. Each step of the workflow can run on a separate machine, enabling parallel execution at scale.

**CloudFormation Sets Up the Environment:**
We use an AWS CloudFormation template to provision the compute environment. This includes setting up the EC2 instances, attaching IAM roles, and configuring queues and resources.

**Tracer Agent is Automatically Installed:**
As part of the CloudFormation provisioning process, each EC2 instance downloads and installs the Tracer agent. This ensures full observability across every job executed by AWS Batch, without requiring any manual setup on the node.

**Session Correlation via Environment Variable:**
To track and group jobs from the same Nextflow pipeline run, we inject the workflow.sessionId as an environment variable into each process. This enables Tracer to associate logs, metrics, and traces across all participating hosts and jobs under one session.

**Fixed Pipeline Name in CloudFormation Template:**
The CloudFormation template used in this setup includes a predefined pipeline name for Tracer. This ensures consistency and simplifies trace grouping across environments.

**Two Testing Options Available:**
You can test this setup either locally (using your own AWS account) or through the Tracer-provided sandbox environment, which is preconfigured for rapid experimentation.


### Option 1: Local Execution
1. Clone or navigate to the [bioinformatics pipelines repository](https://github.com/Tracer-Cloud/tracer-test-pipelines-bioinformatics).
2. Install Nextflow and Spack locally (refer to the repository's README for detailed instructions).
3. Ensure Tracer is installed and configured.
4. Run the following Nextflow command:

```bash
nextflow -c nextflow-config/batch.config run https://github.com/nf-core/rnaseq \
    --outdir output/rnaseq-test \
    -params-file nextflow-config/rnaseq-params.json \
    -profile test
```

This will submit your job to AWS Batch. Pipeline execution should become visible in Grafana.

### Option 2: Sandbox Environment (Currently Only On Tracer AWS Account)

1. Spin up an EC2 Sandbox environment using the provided EC2 launch template.
2. Access the sandbox and switch to the correct user directory:

```bash
sudo su - ubuntu && cd test-bioinformatics-packages
```

3. Initialize Tracer for AWS Batch RNA-seq pipeline:

```bash
tracer init --pipeline-name aws_batch_test \
    --environment sandbox \
    --user-operator vincent \
    --pipeline-type aws_batch_rnaseq
```

4. Run the RNA-seq test pipeline:

```bash
make test_rnaseq_aws_batch
```


## Troubleshooting Local Tests with Cargo

If you encounter the following error while running Cargo tests:

```bash
---- alert stdout ----
thread 'alert' panicked at /path/to/sqlx-postgres/src/testing/mod.rs:78:44:
DATABASE_URL must be set: EnvVar(NotPresent)
```

Set the `DATABASE_URL` environment variable:

```bash
export DATABASE_URL=postgres://tracer-cluster-v2-instance-1.cdgizpzxtdp6.us-east-1.rds.amazonaws.com:5432/tracer_db
```

Then, retry running your tests.

## Running a Local Environment Test

Initialize a local RNA-seq pipeline demo:

```bash
tracer init --pipeline-name rnaseq-demo-123 \
    --environment demo \
    --user-operator vincent \
    --pipeline-type rnaseq
```