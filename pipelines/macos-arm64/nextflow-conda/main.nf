#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "*.fasta"
params.outdir = "results"


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
    MULTIQC_SIM()
    BEDTOOLS_SIM()
}
