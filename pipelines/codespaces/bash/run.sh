#!/bin/bash
# This script has been shared by @ceisenhart to the Tracer team
# Usage: ./bioinformatics_pipeline.sh <fastq1> <fastq2> <output_vcf> [ref_dir]
set -e  # Exit on error

# Validate arguments
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <fastq1> <fastq2> <output_vcf> [ref_dir]"
    exit 1
fi

FASTQ1=$(realpath "$1" 2>/dev/null || { echo "Error: $1 not found"; exit 1; })
FASTQ2=$(realpath "$2" 2>/dev/null || { echo "Error: $2 not found"; exit 1; })
OUTPUT_VCF=$(realpath "$(dirname "$3")/$(basename "$3")" 2>/dev/null || { echo "Error: Invalid output VCF path"; exit 1; })
DATA_DIR=$(dirname "$OUTPUT_VCF")
REF_DIR=${4:-"/Users/ceisenhart/Bioinformatics/data/bwa_indeces"}
REF_FA="hg19_chr19.fa"
BASE_NAME=$(basename "$OUTPUT_VCF" .vcf)

# Validate reference directory and files
if [ ! -d "$REF_DIR" ] || [ ! -f "$REF_DIR/$REF_FA" ]; then
    echo "Error: Reference directory $REF_DIR or $REF_FA not found"
    exit 1
fi

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Pull Docker images
log "Pulling Docker images..."
docker pull pegi3s/fastqc:latest || { echo "Error: Failed to pull pegi3s/fastqc"; exit 1; }
docker pull quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0 || { echo "Error: Failed to pull trim-galore"; exit 1; }
docker pull biocontainers/bwa:v0.7.17_cv1 || { echo "Error: Failed to pull bwa"; exit 1; }
docker pull biocontainers/samtools:v1.9-4-deb_cv1 || { echo "Error: Failed to pull samtools"; exit 1; }
docker pull broadinstitute/gatk:latest || { echo "Error: Failed to pull gatk"; exit 1; }

# 1. FastQC on input FASTQ files
log "Running FastQC on input files..."
docker run --rm -v "$DATA_DIR:/data" pegi3s/fastqc fastqc "/data/$(basename "$FASTQ1")" "/data/$(basename "$FASTQ2")" -o /data || {
    echo "Error: FastQC on input files failed"
    exit 1
}

# 2. TrimGalore
log "Running TrimGalore..."
docker run --rm -v "$DATA_DIR:/data" quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0 trim_galore --paired --dont_gzip "/data/$(basename "$FASTQ1")" "/data/$(basename "$FASTQ2")" -o /data || {
    echo "Error: TrimGalore failed"
    exit 1
}

# 3. FastQC on trimmed files
TRIMMED1="${FASTQ1%.*}_val_1.fq"
TRIMMED2="${FASTQ2%.*}_val_2.fq"
if [ ! -f "$TRIMMED1" ] || [ ! -f "$TRIMMED2" ]; then
    echo "Error: Trimmed files $TRIMMED1 or $TRIMMED2 not found"
    exit 1
fi
log "Running FastQC on trimmed files..."
docker run --rm -v "$DATA_DIR:/data" pegi3s/fastqc fastqc "/data/$(basename "$TRIMMED1")" "/data/$(basename "$TRIMMED2")" -o /data || {
    echo "Error: FastQC on trimmed files failed"
    exit 1
}

# 4. BWA-MEM alignment
log "Running BWA-MEM..."
docker run --rm -v "$DATA_DIR:/data" -v "$REF_DIR:/ref" biocontainers/bwa:v0.7.17_cv1 bwa mem "/ref/$REF_FA" "/data/$(basename "$TRIMMED1")" "/data/$(basename "$TRIMMED2")" -o "/data/$BASE_NAME.sam" || {
    echo "Error: BWA-MEM failed"
    exit 1
}

# 5. Convert SAM to BAM
log "Converting SAM to BAM..."
docker run --rm -v "$DATA_DIR:/data" biocontainers/samtools:v1.9-4-deb_cv1 samtools view -bS "/data/$BASE_NAME.sam" -o "/data/$BASE_NAME.bam" || {
    echo "Error: SAM to BAM conversion failed"
    exit 1
}

# 6. Sort BAM
log "Sorting BAM..."
docker run --rm -v "$DATA_DIR:/data" biocontainers/samtools:v1.9-4-deb_cv1 samtools sort "/data/$BASE_NAME.bam" -o "/data/$BASE_NAME.sorted.bam" || {
    echo "Error: BAM sorting failed"
    exit 1
}

# 7. Index BAM
log "Indexing sorted BAM..."
docker run --rm -v "$DATA_DIR:/data" biocontainers/samtools:v1.9-4-deb_cv1 samtools index "/data/$BASE_NAME.sorted.bam" || {
    echo "Error: BAM indexing failed"
    exit 1
}

# 8. Add read group to BAM
log "Adding read group to BAM..."
docker run --rm -v "$DATA_DIR:/data" biocontainers/samtools:v1.9-4-deb_cv1 samtools addreplacerg -r "ID:RG1\tSM:SAMPLE1\tLB:LIB1\tPL:ILLUMINA" "/data/$BASE_NAME.sorted.bam" -o "/data/$BASE_NAME.rg.bam" || {
    echo "Error: Adding read group failed"
    exit 1
}

# 9. Index read group BAM
log "Indexing read group BAM..."
docker run --rm -v "$DATA_DIR:/data" biocontainers/samtools:v1.9-4-deb_cv1 samtools index "/data/$BASE_NAME.rg.bam" || {
    echo "Error: Read group BAM indexing failed"
    exit 1
}

# 10. Create FASTA index (if not exists)
if [ ! -f "$REF_DIR/$REF_FA.fai" ]; then
    log "Creating FASTA index..."
    docker run --rm -v "$REF_DIR:/ref" biocontainers/samtools:v1.9-4-deb_cv1 samtools faidx "/ref/$REF_FA" || {
        echo "Error: FASTA index creation failed"
        exit 1
    }
fi

# 11. Create FASTA dictionary (if not exists)
if [ ! -f "$REF_DIR/${REF_FA%.fa}.dict" ]; then
    log "Creating FASTA dictionary..."
    docker run --rm -v "$REF_DIR:/ref" broadinstitute/gatk gatk CreateSequenceDictionary -R "/ref/$REF_FA" -O "/ref/${REF_FA%.fa}.dict" || {
        echo "Error: FASTA dictionary creation failed"
        exit 1
    }
fi

# 12. GATK HaplotypeCaller
log "Running GATK HaplotypeCaller..."
docker run --rm -v "$DATA_DIR:/data" -v "$REF_DIR:/ref" broadinstitute/gatk gatk HaplotypeCaller -R "/ref/$REF_FA" -I "/data/$BASE_NAME.rg.bam" -O "/data/$(basename "$OUTPUT_VCF")" || {
    echo "Error: GATK HaplotypeCaller failed"
    exit 1
}

log "Pipeline completed successfully. Output VCF: $OUTPUT_VCF"
