# #!/bin/bash
# set -e

# echo "[*] Updating system packages"
# sudo apt-get update && sudo apt-get install -y \
#     make \
#     curl \
#     wget \
#     unzip \
#     openjdk-11-jdk \
#     git \
#     graphviz \
#     libbz2-dev \
#     liblzma-dev \
#     libcurl4-openssl-dev \
#     libssl-dev \
#     zlib1g-dev \
#     libncurses5-dev \
#     libncursesw5-dev \
#     libreadline-dev

# echo "[*] Installing Nextflow"
# curl -s https://get.nextflow.io | bash
# sudo mv nextflow /usr/local/bin/

# echo "[*] Installing Miniconda and Mamba"
# cd /tmp
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# bash miniconda.sh -b -p $HOME/miniconda
# eval "$($HOME/miniconda/bin/conda shell.bash hook)"
# conda init bash
# source ~/.bashrc
# conda config --set always_yes yes --set changeps1 no
# conda install -c conda-forge mamba

# echo "[*] Creating Conda env with RNASeq tools"
# mamba create -n rnaseq -c bioconda -c conda-forge \
#     fastqc \
#     hisat2 \
#     samtools \
#     subread

# echo 'conda activate rnaseq' >> ~/.bashrc

# echo "[*] Installing Tracer CLI"
# curl -sSL https://install.tracer.cloud/installation-script-development.sh | bash
# source ~/.bashrc

# echo "[✓] Setup complete. Tools installed:"
# which nextflow
# which tracer
# conda info --envs