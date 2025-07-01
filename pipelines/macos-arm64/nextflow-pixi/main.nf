#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "*.fasta"
params.outdir = "results"
params.iterations = 20

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
    
    // Create a range channel for 20 iterations
    iterations = Channel.of(1..params.iterations)
    
    // Run each process 20 times
    STAR_VERSION(iterations)
    MULTIQC_SIM(iterations)
    BEDTOOLS_SIM(iterations)
}
