#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params {
    input = "test_data/*.fasta"
    outdir = "pipelines/linux-arm-ubuntu/results"
}

workflow {
    GET_VERSIONS()
    FASTQC()
}

Channel
    .fromPath(params.input)
    .ifEmpty { error "❌ No input files found at: ${params.input}" }
    .set { fasta_files }

fasta_files.view()

process GET_VERSIONS {
    conda = './environment.yml'

    script:
    """
    fastqc --version
    """
}

process FASTQC {
    conda = './environment.yml'

    input:
        path fasta_file from fasta_files

    output:
        path "*.html"
        path "*.zip"

    script:
    """
    fastqc $fasta_file --outdir ${params.outdir}
    """
}
