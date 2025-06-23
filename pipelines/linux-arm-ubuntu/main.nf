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

    // Execute processes directly — NOT as function calls
    fastqc_process(fasta_files)
    get_versions_process()
}

// FASTQC process
process fastqc_process {
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

// Version checker process
process get_versions_process {
    cpus 1

    script:
    """
    fastqc --version || echo 'fastqc not installed'
    """
}
