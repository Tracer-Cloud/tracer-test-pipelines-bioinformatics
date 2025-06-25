nextflow.enable.dsl = 2

params.outdir = "results"
params.input = "test_data/*.fasta"

include {
    STAR_VERSION,
    SALMON_VERSION,
    FASTA_STATS,
    COUNT_SEQUENCES,
} from '../../shared/nextflow/workflows/fasta-analysis.nf'

workflow {
    // Create input channel for FASTA files
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Run simple version checks for various bioinformatics tools
    fastqc_version()
    star_version()
    samtools_version()
    bwa_version()
    gatk_version()

    // Run shared version checks
    STAR_VERSION()
    SALMON_VERSION()

    // Process FASTA files
    FASTA_STATS(input_ch)
    COUNT_SEQUENCES(input_ch)

    // Collect all version outputs
    fastqc_version.out
        .concat(star_version.out)
        .concat(samtools_version.out)
        .concat(bwa_version.out)
        .concat(gatk_version.out)
        .concat(STAR_VERSION.out)
        .concat(SALMON_VERSION.out)
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

process gatk_version {
    output:
    stdout

    script:
    """
    echo "GATK version:"
    gatk --version || echo "GATK not available"
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
    mkdir -p ${params.outdir}
    cp ${versions_file} tool_versions.txt
    echo "=== Tool Version Summary ===" >> tool_versions.txt
    echo "Pipeline completed at: \$(date)" >> tool_versions.txt
    """
}
