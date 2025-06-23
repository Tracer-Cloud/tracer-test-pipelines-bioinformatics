nextflow.enable.dsl=2

params.input = "test_data/*.fasta"
params.outdir = "results"

workflow {
    fastqc_process(params.input)
    get_versions_process()
}

process fastqc_process {
    tag "$sample"

    input:
    path sample

    output:
    path "${params.outdir}"

    script:
    """
    mkdir -p ${params.outdir}
    fastqc \$sample --outdir ${params.outdir}
    """
}

process get_versions_process {
    output:
    stdout

    script:
    """
    echo "FastQC version:"
    fastqc --version
    """
}
