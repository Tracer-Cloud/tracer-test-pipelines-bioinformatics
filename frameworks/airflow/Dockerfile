FROM apache/airflow:2.7.1

USER root

# 1. Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget bzip2 curl ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 2. Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -u -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set Conda in PATH for all users and shells
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN echo 'export PATH="/opt/conda/bin:$PATH"' > /etc/profile.d/conda.sh


# 3. Set Conda in PATH
ENV PATH="/opt/conda/bin:$PATH"

# 3. Set Conda in PATH
ENV PATH="/opt/conda/bin:$PATH"

# 4. Install mamba FIRST
RUN conda install -y -c conda-forge mamba

# 5. Use conda to configure, then mamba to install
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority flexible && \
    conda install -y python=3.9

# 6. Install all bioinformatics tools with fast resolution
RUN conda install -y -c base \
    fastqc \
    hisat2 \
    samtools \
    bedtools \
    multiqc \
    star \
    rsem \
    stringtie

# 7. Install heavier tools — using bash with .bashrc sourced (optional)
RUN conda install -y -c base \
    star \
    rsem \
    stringtie

# 8. Fix permissions
RUN chgrp -R 0 /opt/conda && chmod -R g+rx /opt/conda/bin

# # 9. Switch to airflow user
USER airflow

# 10. Validate tools
RUN which fastqc && fastqc --version
RUN which STAR && STAR --version
RUN which hisat2 && hisat2 --version
RUN which samtools && samtools --version
RUN which rsem-calculate-expression && rsem-calculate-expression --version
RUN which bedtools && bedtools --version
RUN which multiqc && multiqc --version
