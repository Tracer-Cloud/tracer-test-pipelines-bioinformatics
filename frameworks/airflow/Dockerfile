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

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc && \
    echo "conda activate base" >> /root/.bashrc && \
    /opt/conda/bin/conda init && \
    /opt/conda/bin/conda install -y -c bioconda fastqc

# Ensure conda is activated and in PATH for root during build
RUN bash -c "source /root/.bashrc && conda --version"
RUN bash -c "source /root/.bashrc && which fastqc"
RUN bash -c "source /root/.bashrc && fastqc --version"

# Switch to airflow user for pip installations
USER airflow

# Ensure airflow user has access to conda and fastqc in PATH
ENV PATH="/opt/conda/bin:${PATH}"

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

# Debug: ensure conda works and fastqc is accessible for airflow user
RUN bash -c "source /etc/profile.d/conda.sh && conda --version"
RUN bash -c "which fastqc"
RUN bash -c "fastqc --version"