nextflow.enable.dsl=2

process get_versions_process {
    output:
    stdout into versions_output

    script:
    """
    echo 'fastqc'
    fastqc --version
    """
}

process fastqc_process {
    input:
    path sample_file

    output:
    path "*.html"
    path "*.zip"

    conda:
    'bioconda::fastqc=0.11.9'

    script:
    """
    fastqc $sample_file --outdir ${params.outdir}
    """
}

workflow {
    samples_ch = Channel.fromPath(params.input)
    fastqc_process(samples_ch)
    get_versions_process()
}
