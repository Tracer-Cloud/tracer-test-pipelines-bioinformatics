nextflow.enable.dsl=2

params.input = "test_data/*.fasta"
params.outdir = "pipelines/linux-arm-ubuntu/results"

Channel
    .fromPath(params.input)
    .ifEmpty { error "❌ No input files found at: ${params.input}" }
    .set { fasta_files }

fasta_files.view()

workflow {
    GET_VERSIONS()
    FASTQC(fasta_files)
}

process GET_VERSIONS {
    cpus 1

    script:
    """
    fastqc --version
    """
}

process FASTQC {
    cpus 1

    input:
    path fasta_file

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${fasta_file} --outdir ${params.outdir}
    """
}
