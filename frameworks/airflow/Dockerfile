FROM apache/airflow:2.7.1

USER root

# 1. Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        curl \
        ca-certificates \
        bzip2 \
        libffi-dev \
        libffi7 \
        build-essential \
        software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 2. Install Miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -u -p /opt/conda && \
    rm /tmp/miniconda.sh

# 3. Set environment variables
ENV CONDA_DIR=/opt/conda
ENV PATH="$CONDA_DIR/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"

RUN echo 'export PATH="/opt/conda/bin:$PATH"' > /etc/profile.d/conda.sh

# 4. Install mamba and configure conda channels
RUN conda install -y -c conda-forge mamba && \
    mamba config --add channels bioconda && \
    mamba config --add channels conda-forge && \
    mamba config --set channel_priority strict && \
    mamba install -y python=3.9

# 5. Install bioinformatics tools
RUN mamba install -y \
    fastqc \
    hisat2 \
    samtools \
    bedtools \
    multiqc \
    star \
    rsem \
    stringtie

# 6. Install updated C++ standard lib to fix GLIBCXX issues
RUN mamba install -y libstdcxx-ng && \
    chmod a+rx /opt/conda/lib/libstdc++.so.* && \
    chgrp 0 /opt/conda/lib/libstdc++.so.*

# 7. Ensure bioinformatics tools are executable and shared
RUN chmod a+rx /opt/conda/bin/* && \
    chgrp 0 /opt/conda/bin/*

# 8. Switch back to airflow user
USER airflow
