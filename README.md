<h2 align="left">
Tracer Linux Agent: Observability for Scientific HPC Workloads
</h2>

## Quickstart Tracer

We recommend using the Sandbox Environment for an easy ans quick onboarding experience: https://sandbox.tracer.cloud/

Click the ‘Get started’ button and follow the guided steps—no AWS credentials or setup required.



### 1. Install Tracer With One Line of Code (Already Installed in Codespaces, skip this part)

Install Tracer with this single command:
```bash
curl -sSL https://install.tracer.cloud/installation-script-development.sh | bash && source ~/.bashrc
```
Click the 'Open In Github Codespaces' button to use GitHub Codespaces.

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


### 3. Monitor your Pipeline

Watch your pipeline in action via the Tracer monitoring dashboard, which you access by clicking the ‘Open Grafana Dashboard’ button.
You’ll see real-time execution metrics, stages, and status updates.




<br />



## Table of Contents
- [🔍 Examples](docs/EXAMPLES.md) – Explore real-world use cases 
- [🤝 Contributing](docs/CONTRIBUTING.md) – Join the community and contribute



<br />



## Mission

> *"The goal of Tracer's Rust agent is to equip scientists and engineers with DevOps intelligence to efficiently harness massive computational power for humanity's most critical challenges."*
