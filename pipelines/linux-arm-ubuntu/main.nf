nextflow.enable.dsl = 2

Channel
    .fromPath(params.input)
    .set { fasta_files }

workflow {
    FASTQC(fasta_files)
}

process FASTQC {
    input:
        path fasta_file

    output:
        path "*.html"
        path "*.zip"

    script:
    """
    fastqc $fasta_file
    """
}
