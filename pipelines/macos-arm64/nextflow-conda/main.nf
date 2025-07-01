#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "*.fasta"
params.outdir = "results"
params.iterations = 10

process STAR_VERSION {
    input:
    val iteration
    
    output:
    stdout
    
    script:
    """
    echo "STAR_VERSION iteration ${iteration}"
    STAR --version
    """
}

process MULTIQC_SIM {
    input:
    val iteration
    
    output:
    stdout
    
    script:
    """
    echo "MULTIQC_SIM iteration ${iteration}"
    echo "[multiqc] Simulating long run..."
    multiqc --version
    """
}

process BEDTOOLS_SIM {
    input:
    val iteration
    
    output:
    stdout
    
    script:
    """
    echo "BEDTOOLS_SIM iteration ${iteration}"
    bedtools --version || echo "Bedtools not available"
    """
}

workflow {
    input_ch = Channel.fromPath(params.input)
    
    iterations = Channel.of(1..params.iterations)
    
    STAR_VERSION(iterations)
    MULTIQC_SIM(iterations)
    BEDTOOLS_SIM(iterations)
}
