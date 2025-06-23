nextflow.enable.dsl=2

workflow {
    // Input channel
    Channel
        .fromPath(params.input)
        .ifEmpty { error "❌ No input files found at: ${params.input}" }
        .set { fasta_files }

    // Run processes
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
    fastqc \$fasta_file --outdir ${params.outdir}
    """
}

// Version check process
process get_versions_process {
    script:
    """
    fastqc --version || echo 'fastqc not installed'
    """
}
