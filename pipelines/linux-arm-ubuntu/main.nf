workflow {
    GET_VERSIONS()
    FASTQC()
}

process GET_VERSIONS {
    script:
    """
    fastqc --version
    """
}

// Step 1: Define a channel from input path
Channel
    .fromPath(params.input)
    .ifEmpty { error "❌ No input files found at: ${params.input}" }
    .set { fasta_files }

// Step 2: View matched input files (debugging)
fasta_files.view()

process FASTQC {
    input:
        path fasta_file from fasta_files

    output:
        path "*.html"
        path "*.zip"

    script:
    """
    fastqc $fasta_file
    """
}
