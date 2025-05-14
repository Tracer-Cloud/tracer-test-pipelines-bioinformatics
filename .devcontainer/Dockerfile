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

# Completamento del setup
CMD ["bash"]
