# Nextflow Pipeline on Linux ARM Amazon Linux - Pixi Version

This directory contains a **Pixi-based Nextflow pipeline** for ARM Amazon Linux. The setup and usage are similar to the x86 Ubuntu reference, but adapted for ARM and Amazon Linux.

## Quick Start (Recommended: Pixi)

### Prerequisites

- [Pixi](https://pixi.sh) installed on your system
- ARM Amazon Linux system
- **Minimum 512 MB available memory** (optimized for low-memory ARM instances)

### System Requirements

- **Memory**: 512 MB minimum (pipeline is optimized for low-memory ARM instances)
- **CPU**: 1 core minimum
- **Storage**: ~100 MB for pipeline and dependencies

### Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Run the Pipeline with Pixi

```bash
./run.sh
```

### Memory Optimization

This pipeline has been specifically optimized for ARM instances with limited memory:

- All processes use 512 MB memory maximum
- CPU usage limited to 1 core per process
- Error handling configured to continue on failures
