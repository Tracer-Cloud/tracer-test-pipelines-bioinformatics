#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "*.fasta"
params.outdir = "results"

process STAR_VERSION {
    output:
    stdout
    script:
    """
    STAR --version
    """
}

process MULTIQC_SIM {
    output:
    stdout
    script:
    """
    echo "[multiqc] Simulating long run..."
    multiqc --version
    """
}

process BEDTOOLS_SIM {
    output:
    stdout
    script:
    """
    bedtools --version || echo "Bedtools not available"
    """
}

workflow {
    input_ch = Channel.fromPath(params.input)
    STAR_VERSION()
    MULTIQC_SIM()
    BEDTOOLS_SIM()
}
