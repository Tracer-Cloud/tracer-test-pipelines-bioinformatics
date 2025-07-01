#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = params.outdir ?: "results"
params.iterations = 10

workflow version_check {
    iterations = Channel.of(1..params.iterations)
    
    star_version(iterations)
    multiqc_version(iterations)
    bedtools_version(iterations)
}



process star_version {
    input:
    val iteration
    
    output:
    stdout

    script:
    """
    echo "STAR version (iteration ${iteration}):"
    STAR --version || echo "STAR not available"
    """
}

process multiqc_version {
    input:
    val iteration
    
    output:
    stdout

    script:
    """
    echo "MultiQC version (iteration ${iteration}):"
    multiqc --version || echo "MultiQC not available"
    """
}

process bedtools_version {
    input:
    val iteration
    
    output:
    stdout

    script:
    """
    echo "Bedtools version (iteration ${iteration}):"
    bedtools --version || echo "Bedtools not available"
    """
}
