#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "test_data/*.fasta"
params.outdir = "results"
params.conda_env = "nextflow-minimal"

process STAR_VERSION {
    conda params.conda_env
    
    output:
    stdout
    
    script:
    """
    echo "=== STAR VERSION ==="
    STAR --version
    """
}

process SALMON_VERSION {
    conda params.conda_env
    
    output:
    stdout
    
    script:
    """
    echo "=== SALMON VERSION ==="
    salmon --version
    """
}

process FASTA_STATS {
    publishDir "${params.outdir}/fasta_stats", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    set -euo pipefail
    
    if [ ! -s "${fasta_file}" ]; then
      echo "Input file missing or empty: ${fasta_file}" >&2
      exit 1
    fi
    
    base=\$(basename "${fasta_file}" .fasta)
    seq_count=\$(grep -c '^>' "${fasta_file}")
    total_length=\$(grep -v '^>' "${fasta_file}" | tr -d '\\n' | wc -c)
    
    # Handle division by zero
    if [ "\$seq_count" -eq 0 ]; then
        avg_length=0
    else
        avg_length=\$(echo "scale=2; \$total_length / \$seq_count" | bc -l)
    fi
    
    echo "Sequences: \$seq_count" > "\${base}_stats.txt"
    echo "Total length: \$total_length" >> "\${base}_stats.txt"
    echo "Average length: \$avg_length" >> "\${base}_stats.txt"
    """
}

process COUNT_SEQUENCES {
    publishDir "${params.outdir}/counts", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    set -euo pipefail
    
    if [ ! -s "${fasta_file}" ]; then
      echo "Input file missing or empty: ${fasta_file}" >&2
      exit 1
    fi
    
    base=\$(basename "${fasta_file}" .fasta)
    grep -c '^>' "${fasta_file}" > "\${base}_count.txt"
    """
}

workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Run version checks
    STAR_VERSION()
    SALMON_VERSION()
    
    // Process FASTA files
    FASTA_STATS(input_ch)
    COUNT_SEQUENCES(input_ch)
    
    // Display version information
    STAR_VERSION.out.view { "STAR: $it" }
    SALMON_VERSION.out.view { "Salmon: $it" }
}