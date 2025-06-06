# Minimal Nextflow Pipeline with Conda

A simple, reproducible Nextflow pipeline that runs on Mac using Conda instead of Docker containers.

```bash
# how to run
bash ./setup.sh

nextflow -log logs/stub.log run main.nf -resume
```


# Requirements
### Download and Install Nextflow 
- https://www.nextflow.io/docs/stable/install.html

### Donwload and Install Conda 
- https://docs.conda.io/en/latest/miniconda.html