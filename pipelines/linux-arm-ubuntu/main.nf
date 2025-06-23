// Use Nextflow DSL2
nextflow.enable.dsl=2

// Create channel from input files
Channel
    .fromPath(params.input)
    .ifEmpty { error "❌ No input files found at: ${params.input}" }
    .set { fasta_files }

// Optional: print detected input files for debugging
fasta_files.view()

// Define workflow
workflow {
    GET_VERSIONS()
    FASTQC(fasta_files)
}

// Process to show tool versions
process GET_VERSIONS {
    script:
    """
    fastqc --version
    """
}

// Process to run FastQC on each input file
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
