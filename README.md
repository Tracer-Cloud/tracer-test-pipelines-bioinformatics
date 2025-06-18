# Tracer Linux Agent: Observability for Scientific HPC Workloads

## Quickstart Tracer Sandbox

Follow the instructions on https://sandbox.tracer.cloud/ to get started.

## Pipelines

Below is a list of pipelines available for a quickstart setup in the [pipelines](./pipelines/) directory.

### Install Tracer

Before running any pipeline, you need to install Tracer into your operating system. This is a one-time installation.

For Codespaces, we install the below line of code into the terminal:

```bash
curl -sSL https://install.tracer.cloud/ | bash && source ~/.bashrc
```

### Run tracer init

To launch a new pipeline with a new name, run the following command in the Codespaces terminal:

Make sure the environment has root privileges you're running tracer in, `sudo su` if needed.

```bash
tracer init
```

## Dependencies

Run the dependency install script.

List of all the dependencies:

- Java (OpenJDK 17)
- Python 3
- Miniconda
- Docker
- Nextflow
- Spack

To install all dependencies, simply run:

```bash
bash dependencies_installation.sh
```

This script will install all required tools, including Spack, and set up your environment for running the pipelines.

### Run pipeline

We have pre-installed some pipelines in Codespaces for you to run.

We recommend starting with a simple RNA-seq pipeline in Nextflow:

```bash
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results -resume
```

> ⚠️ This pipeline uses a small dataset for demonstration purposes. Feel free to adapt the dataset or explore other prepared pipelines under the pipelines tab.

### Monitor your Pipeline

Watch your pipeline in action via the Tracer monitoring dashboard, which you access by clicking the 'Open Grafana Dashboard' button.

You'll see real-time execution metrics, stages, and status updates.

## What Is Tracer and Why Use It?

Tracer is a system-level observability platform purpose-built for scientific computing. It combines cutting-edge technological advances with the deep understanding of scientific industries to give insights into their speed and costs. Its one-line install Linux agent and instant dashboards allow for real-time insights into scientific computing environments.

Unlike industry agnostic monitoring agents, Tracer structures DevOps data for scientific pipelines, providing clear visibility into pipeline stages and execution runs. In environments like AWS Batch, where processes and containers are loosely connected, users struggle to understand which processes belong to which pipeline run, and frequently lose logs from failed containers, making debugging difficult.

Tracer solves this by intelligently organizing and labeling pipelines, execution runs, and steps. Because it runs directly on Linux, it requires no code changes and supports any programming language, unlike point solutions that work only with one framework. This makes integration effortless even across multi-workload IT environments, including AlphaFold, Slurm, Airflow, Nextflow and also local Bash scripts.

Architected for regulated industries, it ensures enterprise-grade security, with data never leaving your infrastructure, which is not the case with solutions such as DataDog.

## Key Features

New metrics that help you speed up your pipelines and maximize your budget:

- Time and cost per dataset processed
- Execution duration and bottleneck identification for each pipeline step
- Cost attribution across pipelines, teams, and environments (dev, CI/CD, prod)

Overall, making sense of scientific toolchains with poor/no observability.

## Mission

> "_The goal of Tracer's Rust agent is to equip scientists and engineers with DevOps intelligence to efficiently harness massive computational power for humanity's most critical challenges._"
