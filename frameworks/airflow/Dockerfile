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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# Install Miniconda
RUN apt-get update && apt-get install -y wget bzip2 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "conda activate base" >> /etc/bash.bashrc && \
    /opt/conda/bin/conda init && \
    /opt/conda/bin/conda install -y -c bioconda fastqc && \
    apt-get clean

ENV PATH="/opt/conda/bin:$PATH"


# Architecture detection for system packages
RUN echo "Detecting architecture: $(uname -m)" && \
    if [ "$(uname -m)" = "aarch64" ] || [ "$(uname -m)" = "arm64" ]; then \
        echo "Running on ARM architecture"; \
        # If you need specific ARM system packages, add them here \
        apt-get update && \
        apt-get install -y --no-install-recommends \
            # arm-specific-packages-here \
            procps \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*; \
    else \
        echo "Running on x86_64 architecture"; \
        # If you need specific x86_64 system packages, add them here \
        apt-get update && \
        apt-get install -y --no-install-recommends \
            # x86-specific-packages-here \
            procps \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*; \
    fi

# For bioinformatics tools like fastqc that might need special handling
# Commenting out as these are placeholders - uncomment and modify as needed
# RUN if [ "$(uname -m)" = "aarch64" ] || [ "$(uname -m)" = "arm64" ]; then \
#         echo "Installing ARM-compatible bioinformatics tools"; \
#         # ARM-specific installation method for bioinformatics tools \
#     else \
#         echo "Installing x86_64 bioinformatics tools"; \
#         # Standard installation method for bioinformatics tools \
#         curl -fsSL https://example.com/tool-installer.sh | bash; \
#     fi

# Switch to airflow user for pip installations (required by Airflow)
USER airflow

# Install Python packages as airflow user
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    boto3 \
    pandas \
    matplotlib \
    numpy \
    scipy

# Architecture-specific Python packages if needed
RUN if [ "$(uname -m)" = "aarch64" ] || [ "$(uname -m)" = "arm64" ]; then \
        echo "Installing ARM-compatible Python packages"; \
        # If you have ARM-specific Python package requirements \
        pip install --no-cache-dir \
            # arm-specific-python-packages \
            six; \
    else \
        echo "Installing x86_64 Python packages"; \
        # If you have x86_64-specific Python package requirements \
        pip install --no-cache-dir \
            # x86-specific-python-packages \
            six; \
    fi

# Verify airflow is properly installed
RUN python -c "import airflow; print(f'Airflow version: {airflow.__version__}')"

# Always end with the airflow user for security best practices
USER airflow