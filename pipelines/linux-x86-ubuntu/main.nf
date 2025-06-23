nextflow.enable.dsl = 2

params.input = "test_data/*.fasta"
params.outdir = "results"

workflow {
    Channel
        .fromPath(params.input)
        .ifEmpty { error "❌ No input files found at: ${params.input}" }
        .set { fasta_files }

    fastqc_process(fasta_files)
}

process fastqc_process {
    label 'conda'

    input:
    path fasta_file

    output:
    path "${params.outdir}"

    script:
    """
    mkdir -p ${params.outdir}
    fastqc \$fasta_file --outdir ${params.outdir}
    """
}
