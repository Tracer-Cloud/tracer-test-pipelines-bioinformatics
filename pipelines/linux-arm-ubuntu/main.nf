params.input      = "${params.input}"
params.outdir     = "${params.outdir}"
params.genome_dir = "${params.genome_dir}"

process GET_VERSIONS {
    output:
    stdout into version_log

    script:
    """
    fastqc --version
    STAR --version
    """
}

process FASTQC {
    input:
    path input_files from Channel.fromPath(params.input)
    output:
    path "${input_files.simpleName}_fastqc.html" into qc_files

    script:
    """
    fastqc $input_files
    """
}

process STAR_ALIGN {
    input:
    path input_files from Channel.fromPath(params.input)
    output:
    path "star_output_${input_files.simpleName}.bam"

    script:
    """
    STAR --genomeDir ${params.genome_dir} --readFilesIn $input_files --runThreadN 2 --outFileNamePrefix star_output_
    """
}
