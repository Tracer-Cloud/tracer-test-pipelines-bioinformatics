FROM apache/airflow:2.7.1

USER root
# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        bzip2 \
        curl \
        ca-certificates \
        git \
        procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# # Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "conda activate base" >> /etc/bash.bashrc && \
    /opt/conda/bin/conda init && \
    /opt/conda/bin/conda install -y -c bioconda fastqc

# # Ensure airflow remains accessible after conda setup
# ENV PATH="/home/airflow/.local/bin:/opt/conda/bin:$PATH"

# Architecture detection logging
RUN echo "Detected architecture: $(uname -m)"

# Switch to airflow user for pip installations
USER airflow

# Install Python packages as airflow user
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    boto3 \
    pandas \
    matplotlib \
    numpy \
    scipy \
    six

# Debug: ensure Airflow CLI is accessible
RUN which airflow && airflow version