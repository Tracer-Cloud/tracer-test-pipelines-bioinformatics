nextflow.enable.dsl = 2

process fastqc_version {
    output:
    stdout

    script:
    """
    echo "FastQC version:"
    fastqc --version || echo "FastQC not available"
    """
}

process samtools_version {
    output:
    stdout

    script:
    """
    echo "Samtools version:"
    samtools sort --version || echo "Samtools not available"
    """
}

process star_version {
    output:
    stdout

    script:
    """
    echo "STAR version:"
    STAR --version || echo "STAR not available"
    """
}

process seqtk_version {
    output:
    stdout

    script:
    """
    echo "Seqtk version:"
    seqtk 2>&1 | head -1 || echo "Seqtk not available"
    """
}

workflow {
    fastqc_version()
    samtools_version()
    star_version()
    seqtk_version()
} 