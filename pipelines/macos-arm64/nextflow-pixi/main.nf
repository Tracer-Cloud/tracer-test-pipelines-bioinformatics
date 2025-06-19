#!/usr/bin/env nextflow

// Simple pipeline that processes FASTA files
// Uses FastQC for quality control and basic Unix tools
// Optimized for Pixi environment management

nextflow.enable.dsl = 2

// Parameters
params.input = "*.fasta"
params.outdir = "results"

process STAR_VERSION {
    label 'version_check'

    script:
    """
    echo "=== Tool Versions in Pixi Environment ==="
    echo "FASTQC VERSION:"
    fastqc --version

    echo "STAR VERSION:"
    STAR --version
    
    echo "Environment managed by Pixi"
    """
}

// Process 1: Analyze FASTA files (basic statistics)
process FASTA_STATS {
    publishDir "${params.outdir}/fasta_stats", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_stats.txt"

    script:
    """
    echo "=== FASTA Statistics for ${fasta_file} ===" > ${fasta_file.baseName}_stats.txt
    echo "File: ${fasta_file}" >> ${fasta_file.baseName}_stats.txt
    echo "Generated on: \$(date)" >> ${fasta_file.baseName}_stats.txt
    echo "Environment: Pixi-managed" >> ${fasta_file.baseName}_stats.txt
    echo "" >> ${fasta_file.baseName}_stats.txt

    # Count sequences
    seq_count=\$(grep -c '^>' ${fasta_file})
    echo "Number of sequences: \$seq_count" >> ${fasta_file.baseName}_stats.txt

    # Calculate total length
    total_length=\$(grep -v '^>' ${fasta_file} | tr -d '\\n' | wc -c)
    echo "Total sequence length: \$total_length bp" >> ${fasta_file.baseName}_stats.txt

    # Calculate average length
    if [ \$seq_count -gt 0 ]; then
        avg_length=\$(echo "scale=2; \$total_length / \$seq_count" | bc -l)
        echo "Average sequence length: \$avg_length bp" >> ${fasta_file.baseName}_stats.txt
    fi

    # Show sequence IDs
    echo "" >> ${fasta_file.baseName}_stats.txt
    echo "Sequence IDs:" >> ${fasta_file.baseName}_stats.txt
    grep '^>' ${fasta_file} >> ${fasta_file.baseName}_stats.txt
    """
}

// Process 2: Count sequences in FASTA
process COUNT_SEQUENCES {
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    path fasta_file
    
    output:
    path "${fasta_file.baseName}_count.txt"
    
    script:
    """
    grep -c '^>' ${fasta_file} > ${fasta_file.baseName}_count.txt
    echo "File: ${fasta_file}" >> ${fasta_file.baseName}_count.txt
    echo "Sequences: \$(grep -c '^>' ${fasta_file})" >> ${fasta_file.baseName}_count.txt
    echo "Environment: Pixi-managed" >> ${fasta_file.baseName}_count.txt
    """
}

// Workflow
workflow {
    // Create channel from input files
    input_ch = Channel.fromPath(params.input)
    
    // Log the files being processed
    input_ch.view { "Processing file: $it" }
    
    // Create a copy of the channel for each process
    input_ch_copy = input_ch.collect().flatten()

    STAR_VERSION()
    
    // Run FASTA statistics
    FASTA_STATS(input_ch_copy)

    // Count sequences
    COUNT_SEQUENCES(input_ch_copy)
}
