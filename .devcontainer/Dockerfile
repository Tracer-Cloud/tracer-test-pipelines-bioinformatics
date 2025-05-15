FROM ubuntu:24.04

# Install basic dependencies - cleaned up apt-get commands to reduce image size
RUN apt-get update && apt-get install -y \
  build-essential curl wget git unzip openjdk-17-jdk \
  graphviz make libbz2-dev liblzma-dev libcurl4-openssl-dev \
  libssl-dev zlib1g-dev libncurses5-dev libncursesw5-dev \
  libreadline-dev locales ca-certificates gnupg \
  lsb-release software-properties-common \
  && rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8

# Docker CLI only - rely on host Docker daemon
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" > /etc/apt/sources.list.d/docker.list && \
    apt-get update && \
    apt-get install -y docker-ce-cli && \
    rm -rf /var/lib/apt/lists/*

# Install Tracer
RUN curl -sSL https://install.tracer.cloud/ -o /tmp/install.sh && \
    bash /tmp/install.sh && \
    rm /tmp/install.sh

# Add Tracer to PATH
ENV PATH="/root/.tracer/bin:${PATH}"

# Install Nextflow
RUN wget -O nextflow https://github.com/nextflow-io/nextflow/releases/download/v24.04.2/nextflow && \
    chmod +x nextflow && \
    mv nextflow /usr/local/bin/nextflow


# Install AWS CLI (useful for data transfer to/from S3)
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf aws awscliv2.zip