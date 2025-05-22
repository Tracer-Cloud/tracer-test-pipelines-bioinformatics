#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "data/*.fastq.gz"
params.outdir = "results"

workflow {

    Channel
        .fromPath(params.reads)
        .ifEmpty { error "No input FASTQ files found in: ${params.reads}" }
        .set { ch_reads }

    fastqc(ch_reads)
    trim_galore(ch_reads)
    infer_strand(ch_reads)
    multiqc(Channel.of(params.outdir))
}

process fastqc {
    tag "$sample"
    input:
    path sample

    output:
    path "fastqc/*"

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
    path "*.trimmed.fq.gz"

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

process multiqc {
    input:
    val outdir

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc $outdir -o ./
    """
}