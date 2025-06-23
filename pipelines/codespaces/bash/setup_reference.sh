#!/bin/bash
# Setup script to download reference genome and create directory structure
set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="$SCRIPT_DIR/reference"

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log "Setting up reference genome..."

# Create reference directory
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# Download chromosome 19 from human genome (smaller for testing)
log "Downloading chromosome 19 reference genome..."
if [ ! -f "chr19.fa" ]; then
    # Download chromosome 19 from UCSC
    wget -O chr19.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz"
    gunzip chr19.fa.gz
    log "Downloaded and extracted chr19.fa"
else
    log "chr19.fa already exists, skipping download"
fi

# Create BWA index
log "Creating BWA index..."
if [ ! -f "chr19.fa.bwt" ]; then
    docker run --rm -v "$REF_DIR:/ref" biocontainers/bwa:v0.7.17_cv1 bwa index "/ref/chr19.fa"
    log "Created BWA index"
else
    log "BWA index already exists, skipping"
fi

# Create FASTA index
log "Creating FASTA index..."
if [ ! -f "chr19.fa.fai" ]; then
    docker run --rm -v "$REF_DIR:/ref" biocontainers/samtools:v1.9-4-deb_cv1 samtools faidx "/ref/chr19.fa"
    log "Created FASTA index"
else
    log "FASTA index already exists, skipping"
fi

# Create FASTA dictionary
log "Creating FASTA dictionary..."
if [ ! -f "chr19.dict" ]; then
    docker run --rm -v "$REF_DIR:/ref" broadinstitute/gatk gatk CreateSequenceDictionary -R "/ref/chr19.fa" -O "/ref/chr19.dict"
    log "Created FASTA dictionary"
else
    log "FASTA dictionary already exists, skipping"
fi

log "Reference setup completed successfully!"
log "Reference directory: $REF_DIR"
log "Files created:"
ls -la "$REF_DIR" 