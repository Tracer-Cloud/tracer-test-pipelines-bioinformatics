# Tracer Linux Agent: Observability for Scientific HPC Workloads

## Quickstart Tracer Sandbox

Get started instantly by visiting [sandbox.tracer.cloud](https://sandbox.tracer.cloud/).

## Getting Started

### 1. Install Tracer

Install Tracer on your operating system (one-time installation):

```bash
curl -sSL https://install.tracer.cloud/ | bash && source ~/.bashrc
```

### 2. Initialize Tracer

Launch Tracer by running: (Root privileges required)

```bash
sudo tracer init
```

### 3. Choose a Pipeline

We provide several ready-to-use pipelines for different environments in the [pipelines](./pipelines/) directory. Navigate to the pipelines directory and follow the README instructions to choose and run a pipeline suitable for your environment.

## Monitor your Pipeline

Track your pipeline's progress through the Tracer monitoring dashboard, accessible via the 'Open Grafana Dashboard' button.

The dashboard provides real-time insights into:

- Execution metrics
- Pipeline stages
- Status updates

## What Is Tracer and Why Use It?

Tracer is a cutting-edge system-level observability platform specifically designed for scientific computing. It combines advanced technology with deep industry knowledge to provide comprehensive insights into performance and costs. With its simple one-line Linux agent installation and intuitive dashboards, Tracer delivers immediate visibility into scientific computing environments.

Unlike general-purpose monitoring tools, Tracer is purpose-built for scientific pipelines, offering clear visibility into pipeline stages and execution runs. This is particularly valuable in environments like AWS Batch, where tracking processes across containers can be challenging and failed container logs are often lost.

Tracer excels by:

- Intelligently organizing and labeling pipelines, execution runs, and steps
- Running directly on Linux without requiring code modifications
- Supporting any programming language
- Enabling seamless integration across diverse IT environments (AlphaFold, Slurm, Airflow, Nextflow, and local Bash scripts)

Built with enterprise security in mind, Tracer ensures your data never leaves your infrastructure - a key advantage over solutions like DataDog.

## Key Features

Optimize your pipelines with powerful metrics:

- Time and cost per dataset processed
- Execution duration and bottleneck identification for each pipeline step
- Cost attribution across pipelines, teams, and environments (dev, CI/CD, prod)

These insights help make sense of complex scientific toolchains that traditionally lack proper observability.

## Mission

> "_The goal of Tracer's Rust agent is to equip scientists and engineers with DevOps intelligence to efficiently harness massive computational power for humanity's most critical challenges._"
