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

    // Run FastQC process on each FASTA file
    FASTQC_PROCESS(fasta_files)

    // Print FastQC version
    GET_VERSIONS_PROCESS()
}

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

process GET_VERSIONS_PROCESS {
    echo true
    cpus 1

    script:
    """
    fastqc --version || echo 'fastqc not installed'
    """
}
