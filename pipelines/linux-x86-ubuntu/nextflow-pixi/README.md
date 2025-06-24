# Linux x86-64 Ubuntu Nextflow Pipeline with Pixi

This directory contains a Nextflow pipeline setup for Linux x86-64 Ubuntu systems using Pixi package management.

## Quick Start

1. **Run the setup script:**
   ```bash
   ./run.sh
   ```

2. **Activate conda in your current shell session:**
   ```bash
   source ~/.bashrc
   ```

   Or alternatively, use the helper script:
   ```bash
   source activate-conda.sh
   ```

3. **Verify conda is available:**
   ```bash
   conda --version
   ```

## Files

- `run.sh` - Main setup script that installs Java, Miniconda, and Nextflow
- `activate-conda.sh` - Helper script to activate conda in current shell session
- `main.nf` - Nextflow pipeline definition
- `nextflow.config` - Nextflow configuration
- `test_data/` - Sample test data for the pipeline

## Troubleshooting

### Conda command not found after running run.sh

This is expected behavior. The `run.sh` script installs conda and configures it, but you need to activate it in your shell session:

**Option 1:** Source your bashrc file
```bash
source ~/.bashrc
```

**Option 2:** Use the activation helper script
```bash
source activate-conda.sh
```

**Option 3:** Restart your terminal/shell session

### Why does this happen?

The setup script modifies your `~/.bashrc` file to include conda in your PATH, but these changes only take effect in new shell sessions or when you explicitly source the bashrc file.

## Usage

Once conda is activated, you can use it normally:
```bash
conda --version
conda list
conda create -n myenv python=3.9
conda activate myenv
```

### Running nf-core pipelines

After setting up the environment, you can run nf-core pipelines. For example, to run the RNA-seq pipeline:

```bash
# Make sure you're in the nextflow-pixi directory
cd pipelines/linux-x86-ubuntu/nextflow-pixi

# Run the nf-core rnaseq pipeline with test data
nextflow run nf-core/rnaseq -c custom.config -profile docker,test --outdir results -resume
```

The `custom.config` file is included in this directory and contains:
- `process.errorStrategy = 'ignore'` - Allows the pipeline to continue even if some processes fail