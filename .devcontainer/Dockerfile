FROM ubuntu:24.04

# Installare le dipendenze di base
RUN apt-get update && apt-get install -y \
  build-essential \
  curl \
  wget \
  git \
  unzip \
  openjdk-11-jdk \
  graphviz \
  make \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  zlib1g-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libreadline-dev \
  locales \
  ca-certificates

# Genera la locale per evitare warning
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8

# Installazione di Nextflow
RUN curl -s https://get.nextflow.io | bash \
  && mv nextflow /usr/local/bin/

# Installazione di Miniconda (Conda) e Mamba
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
  && bash miniconda.sh -b -p $HOME/miniconda \
  && eval "$($HOME/miniconda/bin/conda shell.bash hook)" \
  && conda init bash \
  && conda config --set always_yes yes --set changeps1 no \
  && conda install -c conda-forge mamba

# Crea l'ambiente Conda con i tool RNASeq
RUN mamba create -n rnaseq -c bioconda -c conda-forge \
  fastqc \
  hisat2 \
  samtools \
  subread \
  && echo 'conda activate rnaseq' >> ~/.bashrc

# Installazione di Tracer (senza inizializzarlo)
RUN curl -sSL https://install.tracer.cloud/installation-script-development.sh | bash

# Imposta il PATH di conda, nextflow, e tracer per essere sicuri che siano disponibili nel container
ENV PATH="/root/miniconda/bin:$PATH"
ENV PATH="/usr/local/bin:$PATH"

# Rendi l'ambiente Conda attivo
RUN echo "conda activate rnaseq" >> ~/.bashrc

# Completamento del setup
CMD ["bash"]
