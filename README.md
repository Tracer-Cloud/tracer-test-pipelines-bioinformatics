<h1 align="left">
Tracer Linux Agent: Observability for Scientific HPC Workloads
</h1>

## Quickstart Tracer

We recommend using the Sandbox Environment for an easy ans quick onboarding experience: https://sandbox.tracer.cloud/

Click the ‘Get started’ button and follow the guided steps. You don't need a local instance and can use our web instance. 



### 1. Install Tracer With One Line of Code (Already Installed in Codespaces, skip this part)

Install Tracer with this single command:
```bash
curl -sSL https://install.tracer.cloud/installation-script-development.sh | bash && source ~/.bashrc
```
Click the 'Open In Github Codespaces' button to use GitHub Codespaces.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/Tracer-Cloud/tracer-test-pipelines-bioinformatics)

Once in Codespaces, the environment comes with:
Tracer pre-installed and Docker running a minimal Nextflow example. Here, you need to run the tracer init command showcased in the next step.



### 2. Initialize a Pipeline

Set up your RNA-seq pipeline by running the following command and run Tracer:
```bash
tracer init --pipeline-name demo_username --environment demo --pipeline-type rnaseq --user-operator user_email --is-dev false 
 ```
Then you need to run a Nextflow command example.

### 3. Run a simple nextflow RNASeq Pipeline
You can run a simple RNASeq pipeline by using this command
```bash
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results -resume
```

> ⚠️ **Warning:** This is a small pipeline sample, just for demo purposes, so you will not see many tools, to understand the full potential of Tracer, you can install it in your Ubuntu machine, and run your favourite pipeline.


### 4. Monitor your Pipeline

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




## Table of Contents
- [🔍 Examples](docs/EXAMPLES.md) – Explore real-world use cases 
- [🤝 Contributing](docs/CONTRIBUTING.md) – Join the community and contribute



<br />



## Mission

> *"The goal of Tracer's Rust agent is to equip scientists and engineers with DevOps intelligence to efficiently harness massive computational power for humanity's most critical challenges."*
