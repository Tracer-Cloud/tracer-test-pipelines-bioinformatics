<h1 align="left">
Tracer Linux Agent: Observability for Scientific HPC Workloads
</h1>

## Quickstart Tracer Sandbox

Follow the instructions on https://sandbox.tracer.cloud/ to get started.

### Install Tracer

Before running any pipeline, you need to install Tracer into your operating system. This is a one-time installation.

For Codespaces, we install the below line of code into the terminal:

```bash
curl -sSL https://install.tracer.cloud/ | bash && source ~/.bashrc
```

### Run tracer init

To launch a new pipeline with a new name, run the following command in the Codespaces terminal:

```bash
tracer init
```

### Run pipeline

We have pre-installed some pipelines in the Codespaces for you to run.
We would recommend to start with a simple rnaseq pipeline in Nextflow:

```bash
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results -resume
```

> ⚠️ This pipeline uses a small dataset for demo purposes. Feel free to adapt the dataset or explore other prepard pipelines under the pipelines tab.

### Monitor your Pipeline

Watch your pipeline in action via the Tracer monitoring dashboard, which you access by clicking the ‘Open Grafana Dashboard’ button.
You’ll see real-time execution metrics, stages, and status updates.

<br />

## What Is Tracer and Why Use It?

- Tracer is a system-level observability platform purpose-built for scientific computing. It combines cutting-edge technological advances withthe deep understanding of scientific industries to give insights into their speed and costs.
  Its one-line install Linux agent and instant dashboards allow for real-time insights into scientific computing environments.

- Unlike industry agnostic monitoring agents, Tracer structures DevOps data for scientific pipelines, providing clear visibility into pipeline stages and execution runs. In environments like AWS Batch, where processes and containers are loosely connected, users struggle to understand which processes belong to which pipeline run, and frequently lose logs from failed containers, making debugging difficult.

- Tracer solves this by intelligently organizing and labeling pipelines, execution runs, and steps. Because it runs directly on Linux, it requires no code changes and supports any programming language, unlike point solutions that work only with one framework. This makes integration effortless even across multi-workload IT environments, including AlphaFold, Slurm, Airflow, Nextflow and also local Bash scripts.

- Architected for regulated industries, it ensures enterprise-grade security, with data never leaving your infrastructure, which is not the case with solutions such as DataDog.

<br />

## Key Features

New metrics that help you speed up your pipelines and maximize your budget:

- Time and cost per dataset processed
- Execution duration and bottleneck identification for each pipeline step
- Cost attribution across pipelines, teams, and environments (dev, CI/CD, prod)
  Overall, making sense of scientific toolchains with poor/no observability.

<br />

## Mission

> _"The goal of Tracer's Rust agent is to equip scientists and engineers with DevOps intelligence to efficiently harness massive computational power for humanity's most critical challenges."_
