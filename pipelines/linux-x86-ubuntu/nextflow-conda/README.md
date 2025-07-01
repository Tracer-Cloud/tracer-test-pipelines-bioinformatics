# Nextflow Pipeline on Linux x86 using Conda

# Linux x86-64 Ubuntu Nextflow Pipeline - Conda Version

This directory contains a **Conda-based Nextflow pipeline** for traditional package management. Uses bash scripts and conda environments for tool installation.

## 🚀 Quick Start - Conda Pipeline

### Prerequisites

- Linux x86-64 system
- sudo access for package installation

### One-Command Setup

1. **Run the setup script:**

   ```bash
   ./run.sh
   ```

2. **Activate conda in your current shell session:**

   ```bash
   source ~/.bashrc
   ```

3. **Verify conda is available:**

   ```bash
   conda activate linux-x86-ubuntu-minimal
   ```

4. ** Run the script manually**

```bash
nextflow -log logs/nextflow.log run main.nf --outdir results
```
