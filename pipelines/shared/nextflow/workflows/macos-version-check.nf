#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = params.outdir ?: "results"
params.iterations = 20

workflow version_check {
    // Create channel with 5 iterations
    iterations = Channel.of(1..params.iterations)
    
    fastqc_version(iterations)
    star_version(iterations)
    multiqc_version(iterations)
    bedtools_version(iterations)

    fastqc_version.out
        .concat(star_version.out)
        .concat(multiqc_version.out)
        .concat(bedtools_version.out)
        .collectFile(name: 'tool_versions.txt', newLine: true)
        .set { all_versions }
}

process fastqc_version {
    input:
    val iteration
    
    output:
    stdout

    script:
    """
    echo "FastQC version (iteration ${iteration}):"
    fastqc --version || echo "FastQC not available"
    """
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
