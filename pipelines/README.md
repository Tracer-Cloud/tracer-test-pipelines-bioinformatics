# Tracer Bioinformatics Test Pipelines

This repository contains validated test pipelines for various bioinformatics platforms and environments, organized by architecture and package management solutions.

## Supported Platforms

Choose your platform below to find validated examples with different package managers:

- [macOS ARM64 (Apple Silicon M1/M2)](macos-arm64/)
- [macOS Intel (x86_64)](macos-intel-x86/)
- [Linux ARM (Ubuntu/CentOS)](linux-arm-ubuntu/)
- [Linux x86 (Ubuntu/CentOS)](linux-x86-ubuntu/)

Navigate to your platform's directory to find working examples with supported package managers.

## Package Manager Installation

We support the following package managers. Choose one based on your needs:

### Pixi (Recommended)

Install Pixi using the official installer:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

For detailed installation options and troubleshooting, visit [Pixi Installation Guide](https://pixi.sh/latest/installation/).

### Conda

Choose one of the following installers:

- **Miniconda**: Minimal installer (~400MB) - Recommended for most users
- **Miniforge**: Community-maintained, conda-forge focused

Installation guides:

- [macOS Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)
- [Linux Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

### Spack

For Spack installation and setup, visit [Spack Tutorial](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html).

## Getting Started

1. Click on your platform's link above
2. Choose an example with your preferred package manager
3. Follow the README instructions in that directory to run the example
