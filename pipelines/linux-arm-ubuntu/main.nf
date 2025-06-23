nextflow.enable.dsl = 2

workflow {
    GET_VERSIONS()
    FASTQC()
    STAR_ALIGN()
}

process GET_VERSIONS {
    script:
    """
    fastqc --version
    """
}

process FASTQC {
    input:
        path fastq_files from file(params.input)

    output:
        path "*.html"
        path "*.zip"

    script:
    """
    fastqc $fastq_files
    """
}

process STAR_ALIGN {
    input:
        path fastq_files from file(params.input)

    output:
        path "*.bam"

    script:
    """
    STAR --genomeDir ${params.genome_dir} --readFilesIn $fastq_files --runThreadN ${task.cpus}
    """
}
