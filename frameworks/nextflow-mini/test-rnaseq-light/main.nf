#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "data/*.fastq.gz"

workflow {
    Channel
        .fromPath(params.reads)
        .ifEmpty { error "No input FASTQ files found in: ${params.reads}" }
        .set { ch_reads }

    fastqc(ch_reads)
    trim_galore(ch_reads)
    infer_strand(ch_reads)
}

process fastqc {
    tag "$sample"
    input:
    path sample

    output:
    path "fastqc"

    script:
    """
    mkdir -p fastqc
    fastqc $sample --outdir fastqc
    """
}

process trim_galore {
    tag "$sample"
    input:
    path sample

    output:
    path "*_trimmed.fq.gz"

    script:
    """
    trim_galore --cores 1 --gzip $sample
    """
}

process infer_strand {
    tag "$sample"
    input:
    path sample

    output:
    path "strandness.txt"

    script:
    """
    infer_experiment.py -i $sample -r dummy.bed > strandness.txt
    """
}