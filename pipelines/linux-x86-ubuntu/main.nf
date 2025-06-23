nextflow.enable.dsl=2

params.input = "$params.input"
params.outdir = "$params.outdir"

workflow {
    Channel
        .fromPath(params.input)
        .ifEmpty { error "❌ No input files found at: ${params.input}" }
        .set { fasta_files }

    fastqc_process(fasta_files)
    get_versions_process()
}

process fastqc_process {
    label 'conda'
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

process get_versions_process {
    label 'conda'
    cpus 1

    script:
    """
    fastqc --version || echo 'fastqc not installed'
    """
}
