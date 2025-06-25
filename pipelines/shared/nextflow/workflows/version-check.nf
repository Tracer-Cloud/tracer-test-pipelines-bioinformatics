nextflow.enable.dsl = 2

params.outdir = "results"

workflow {
    // Run simple version checks for various bioinformatics tools
    fastqc_version()
    star_version()
    samtools_version()
    bwa_version()

    // Collect all version outputs
    fastqc_version.out
        .concat(star_version.out)
        .concat(samtools_version.out)
        .concat(bwa_version.out)
        .collectFile(name: 'tool_versions.txt', newLine: true)
        .set { all_versions }

    // Save results
    save_results(all_versions)
}

process fastqc_version {
    output:
    stdout

    script:
    """
    echo "FastQC version:"
    fastqc --version || echo "FastQC not available"
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

process samtools_version {
    output:
    stdout

    script:
    """
    echo "Samtools version:"
    samtools --version || echo "Samtools not available"
    """
}

process bwa_version {
    output:
    stdout

    script:
    """
    echo "BWA version:"
    bwa 2>&1 | head -3 || echo "BWA not available"
    """
}

process save_results {
    publishDir params.outdir, mode: 'copy'

    input:
    path versions_file

    output:
    path "tool_versions.txt"

    script:
    """
    cp $versions_file tool_versions.txt.tmp
    echo "=== Tool Version Summary ===" >> tool_versions.txt.tmp
    echo "Pipeline completed at: $(date)" >> tool_versions.txt.tmp
    mv tool_versions.txt.tmp tool_versions.txt
    """
} 