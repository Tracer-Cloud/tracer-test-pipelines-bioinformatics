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

process GET_VERSIONS {
    echo true
    cpus 1

    """
    fastqc --version || echo 'fastqc not installed'
    """
}

process FASTQC {
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
