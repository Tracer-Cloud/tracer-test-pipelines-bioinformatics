nextflow.enable.dsl=2

workflow {
    main:
        fastqc_process(params.input)
        get_versions_process()

    emit:
        fastqc_process.out
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
    fastqc $sample --outdir ${params.outdir}
    """
}

process get_versions_process {
    output:
    stdout into version_output

    script:
    """
    echo "FastQC version:"
    fastqc --version
    """
}
