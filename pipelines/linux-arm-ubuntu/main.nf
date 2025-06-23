nextflow.enable.dsl=2

params {
    input = "test_data/*.fasta"
    outdir = "pipelines/linux-arm-ubuntu/results"
}

workflow {

    Channel
        .fromPath(params.input)
        .ifEmpty { error "❌ No input files found at: ${params.input}" }
        .set { fasta_files }

    GET_VERSIONS()
    FASTQC(fasta_files)
}

// Define GET_VERSIONS as a workflow component
workflow GET_VERSIONS {
    main:
    GET_VERSIONS_PROCESS()
}

// Define FASTQC as a workflow component
workflow FASTQC {
    take:
    fasta_files

    main:
    FASTQC_PROCESS(fasta_files)
}

// Process to get FastQC version
process GET_VERSIONS_PROCESS {
    echo true
    cpus 1

    """
    fastqc --version || echo 'fastqc not installed'
    """
}

// FASTQC process
process FASTQC_PROCESS {
    publishDir params.outdir, mode: 'copy'

    input:
    path fasta_file

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc $fasta_file --outdir ${params.outdir}
    """
}
