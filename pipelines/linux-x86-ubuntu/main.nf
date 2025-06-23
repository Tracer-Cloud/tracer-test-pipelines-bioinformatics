nextflow.enable.dsl = 2

workflow {
    get_versions_process()
    fastqc_process()
}

process get_versions_process {
    label 'conda'

    output:
    stdout into versions_out

    script:
    """
    fastqc --version
    """
}

process fastqc_process {
    label 'conda'

    input:
    path sample from file(params.input)

    output:
    path "${params.outdir}"

    script:
    """
    mkdir -p ${params.outdir}
    fastqc \$sample --outdir ${params.outdir}
    """
}
