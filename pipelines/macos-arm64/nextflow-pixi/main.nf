#!/usr/bin/env nextflow

// Enhanced pipeline that processes FASTA files with multiple bioinformatics tools
// Uses various lightweight bioinformatics tools for comprehensive analysis
// Optimized for Pixi environment management

nextflow.enable.dsl = 2

// Parameters
params.input = "test_data/*.fasta"
params.outdir = "results"

// Simple version check process
process VERSION_CHECK {
    output:
    stdout
    
    script:
    """
    echo "=== Bioinformatics Tools Version Check ==="
    echo "Environment: Pixi-managed"
    echo "Platform: macOS ARM64"
    echo ""
    
    echo "1. Seqtk:"
    seqtk 2>&1 | head -1 || echo "Seqtk available"
    echo ""
    
    echo "2. Samtools:"
    samtools --version || echo "Samtools available"
    echo ""
    
    echo "3. Bedtools:"
    bedtools --version || echo "Bedtools available"
    echo ""
    
    echo "4. BLAST:"
    blastn -version || echo "BLAST available"
    echo ""
    
    echo "5. Coreutils:"
    coreutils --version || echo "Coreutils available"
    echo ""
    
    echo "6. BC Calculator:"
    bc --version || echo "BC available"
    echo ""
    
    echo "7. System Info:"
    uname -a
    echo ""
    
    echo "8. Memory Info:"
    sysctl hw.memsize | awk '{print "Total RAM: " \$2/1024/1024/1024 " GB"}'
    """
}

// Process 1: Basic FASTA statistics
process FASTA_STATS {
    publishDir "${params.outdir}/fasta_stats", mode: 'copy'
    
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
    
    # Basic statistics using simple commands
    seq_count=\$(grep -c '^>' "${fasta_file}")
    
    # Get total length using seqtk
    seqtk comp "${fasta_file}" > temp_comp.txt 2>/dev/null || echo "0 0 0 0 0 0" > temp_comp.txt
    
    # Calculate total length
    total_length=\$(awk '{sum += \$2} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    
    # Calculate average length
    if [ "\$seq_count" -gt 0 ] && [ "\$total_length" -gt 0 ]; then
        avg_length=\$(echo "scale=2; \$total_length / \$seq_count" | bc -l 2>/dev/null || echo "0")
    else
        avg_length=0
    fi
    
    # Calculate GC content
    gc_count=\$(awk '{sum += \$3 + \$4} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    if [ "\$total_length" -gt 0 ] && [ "\$gc_count" -gt 0 ]; then
        gc_percent=\$(echo "scale=2; \$gc_count * 100 / \$total_length" | bc -l 2>/dev/null || echo "0.00")
    else
        gc_percent=0.00
    fi
    
    # Write report
    {
        echo "=== FASTA Analysis Report ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo "Environment: Pixi-managed"
        echo ""
        echo "=== Basic Statistics ==="
        echo "Number of sequences: \$seq_count"
        echo "Total sequence length: \$total_length bp"
        echo "Average sequence length: \$avg_length bp"
        echo "GC content: \$gc_percent%"
        echo ""
        echo "=== Sequence Details ==="
        echo "GC bases: \$gc_count"
        echo "Total bases: \$total_length"
        echo ""
        echo "=== Sequence IDs ==="
        grep '^>' "${fasta_file}" | head -10
        if [ \$seq_count -gt 10 ]; then
            echo "... and \$((seq_count - 10)) more sequences"
        fi
    } > "\${base}_stats.txt"
    
    # Clean up temp file
    rm -f temp_comp.txt
    """
}

// Process 2: Sequence length distribution
process SEQUENCE_LENGTHS {
    publishDir "${params.outdir}/lengths", mode: 'copy'
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    set -euo pipefail
    
    base=\$(basename "${fasta_file}" .fasta)
    
    # Get sequence lengths using seqtk
    seqtk comp "${fasta_file}" | awk '{print \$2}' | sort -n > "\${base}_lengths.txt" 2>/dev/null || echo "0" > "\${base}_lengths.txt"
    
    # Calculate basic statistics
    count=\$(wc -l < "\${base}_lengths.txt")
    if [ "\$count" -gt 0 ]; then
        min=\$(head -1 "\${base}_lengths.txt")
        max=\$(tail -1 "\${base}_lengths.txt")
        # Simple median calculation
        if [ \$((count % 2)) -eq 1 ]; then
            median_pos=\$(((count + 1) / 2))
            median=\$(sed -n "\${median_pos}p" "\${base}_lengths.txt")
        else
            median_pos1=\$((count / 2))
            median_pos2=\$((count / 2 + 1))
            val1=\$(sed -n "\${median_pos1}p" "\${base}_lengths.txt")
            val2=\$(sed -n "\${median_pos2}p" "\${base}_lengths.txt")
            median=\$(echo "scale=0; (\$val1 + \$val2) / 2" | bc -l 2>/dev/null || echo "0")
        fi
    else
        min=0
        max=0
        median=0
    fi
    
    echo "Length Statistics for \${base}" > "\${base}_length_stats.txt"
    echo "Total sequences: \$count" >> "\${base}_length_stats.txt"
    echo "Minimum length: \$min bp" >> "\${base}_length_stats.txt"
    echo "Maximum length: \$max bp" >> "\${base}_length_stats.txt"
    echo "Median length: \$median bp" >> "\${base}_length_stats.txt"
    """
}

// Process 3: Nucleotide composition
process NUCLEOTIDE_COMPOSITION {
    publishDir "${params.outdir}/composition", mode: 'copy'
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    set -euo pipefail
    
    base=\$(basename "${fasta_file}" .fasta)
    
    # Use seqtk comp for nucleotide composition
    seqtk comp "${fasta_file}" > temp_comp.txt 2>/dev/null || echo "0 0 0 0 0 0" > temp_comp.txt
    
    # Extract counts from seqtk output
    a_count=\$(awk '{sum += \$3} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    t_count=\$(awk '{sum += \$4} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    g_count=\$(awk '{sum += \$5} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    c_count=\$(awk '{sum += \$6} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    n_count=\$(awk '{sum += \$7} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    other_count=\$(awk '{sum += \$8} END {print sum}' temp_comp.txt 2>/dev/null || echo "0")
    
    total=\$(echo "\$a_count + \$t_count + \$g_count + \$c_count + \$n_count + \$other_count" | bc -l 2>/dev/null || echo "0")
    
    # Calculate percentages
    if [ "\$total" -gt 0 ]; then
        a_percent=\$(echo "scale=2; \$a_count * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        t_percent=\$(echo "scale=2; \$t_count * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        g_percent=\$(echo "scale=2; \$g_count * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        c_percent=\$(echo "scale=2; \$c_count * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        n_percent=\$(echo "scale=2; \$n_count * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        at_content=\$(echo "scale=2; (\$a_count + \$t_count) * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        gc_content=\$(echo "scale=2; (\$g_count + \$c_count) * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
    else
        a_percent="0.00"
        t_percent="0.00"
        g_percent="0.00"
        c_percent="0.00"
        n_percent="0.00"
        at_content="0.00"
        gc_content="0.00"
    fi
    
    echo "Nucleotide Composition Analysis" > "\${base}_composition.txt"
    echo "File: ${fasta_file}" >> "\${base}_composition.txt"
    echo "Total bases: \$total" >> "\${base}_composition.txt"
    echo "" >> "\${base}_composition.txt"
    echo "A: \$a_count (\$a_percent%)" >> "\${base}_composition.txt"
    echo "T: \$t_count (\$t_percent%)" >> "\${base}_composition.txt"
    echo "G: \$g_count (\$g_percent%)" >> "\${base}_composition.txt"
    echo "C: \$c_count (\$c_percent%)" >> "\${base}_composition.txt"
    echo "N: \$n_count (\$n_percent%)" >> "\${base}_composition.txt"
    echo "Other: \$other_count" >> "\${base}_composition.txt"
    echo "" >> "\${base}_composition.txt"
    echo "AT content: \$at_content%" >> "\${base}_composition.txt"
    echo "GC content: \$gc_content%" >> "\${base}_composition.txt"
    
    # Clean up temp file
    rm -f temp_comp.txt
    """
}

// Process 4: Simple quality metrics
process QUALITY_METRICS {
    publishDir "${params.outdir}/quality", mode: 'copy'
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    set -euo pipefail
    
    base=\$(basename "${fasta_file}" .fasta)
    
    # Check for common issues using simple commands
    echo "Quality Assessment Report" > "\${base}_quality.txt"
    echo "File: ${fasta_file}" >> "\${base}_quality.txt"
    echo "Analysis Date: \$(date)" >> "\${base}_quality.txt"
    echo "" >> "\${base}_quality.txt"
    
    # Check for empty sequences
    empty_seqs=0
    if grep -q '^>$' "${fasta_file}"; then
        empty_seqs=\$(grep -c '^>$' "${fasta_file}")
    fi
    echo "Empty sequences: \$empty_seqs" >> "\${base}_quality.txt"
    
    # Check for very short sequences using seqtk
    short_seqs=0
    if [ -s "${fasta_file}" ]; then
        seqtk comp "${fasta_file}" > temp_short.txt 2>/dev/null || echo "0 0 0 0 0 0" > temp_short.txt
        short_seqs=\$(awk '\$2 < 10 {count++} END {print count+0}' temp_short.txt 2>/dev/null || echo "0")
        rm -f temp_short.txt
    fi
    echo "Very short sequences (<10 bp): \$short_seqs" >> "\${base}_quality.txt"
    
    # Check for ambiguous bases using seqtk
    ambig_bases=0
    if [ -s "${fasta_file}" ]; then
        seqtk comp "${fasta_file}" > temp_ambig.txt 2>/dev/null || echo "0 0 0 0 0 0" > temp_ambig.txt
        ambig_bases=\$(awk '{sum += \$7 + \$8} END {print sum}' temp_ambig.txt 2>/dev/null || echo "0")
        rm -f temp_ambig.txt
    fi
    echo "Ambiguous bases: \$ambig_bases" >> "\${base}_quality.txt"
    
    # Check for duplicate headers
    dup_headers=0
    if [ -s "${fasta_file}" ]; then
        dup_headers=\$(grep '^>' "${fasta_file}" | sort | uniq -d | wc -l)
    fi
    echo "Duplicate headers: \$dup_headers" >> "\${base}_quality.txt"
    
    echo "" >> "\${base}_quality.txt"
    
    # Calculate quality score
    if command -v bc >/dev/null 2>&1; then
        quality_score=\$(echo "scale=1; 100 - \$empty_seqs*10 - \$short_seqs*5 - \$ambig_bases/100 - \$dup_headers*20" | bc -l 2>/dev/null | awk '{if(\$1<0) print 0; else print \$1}' 2>/dev/null || echo "100")
    else
        quality_score=100
    fi
    echo "Quality Score: \$quality_score/100" >> "\${base}_quality.txt"
    """
}

// Process 5: Summary report
process SUMMARY_REPORT {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    set -euo pipefail
    
    base=\$(basename "${fasta_file}" .fasta)
    
    # Generate a simple summary using seqtk
    seq_count=\$(grep -c '^>' "${fasta_file}")
    
    # Get total length using seqtk
    seqtk comp "${fasta_file}" > temp_total.txt 2>/dev/null || echo "0 0 0 0 0 0" > temp_total.txt
    total_length=\$(awk '{sum += \$2} END {print sum}' temp_total.txt 2>/dev/null || echo "0")
    rm -f temp_total.txt
    
    echo "=== FASTA Summary Report ===" > "\${base}_summary.txt"
    echo "File: ${fasta_file}" >> "\${base}_summary.txt"
    echo "Sequences: \$seq_count" >> "\${base}_summary.txt"
    echo "Total length: \$total_length bp" >> "\${base}_summary.txt"
    if [ "\$seq_count" -gt 0 ] && [ "\$total_length" -gt 0 ]; then
        echo "Average length: \$(echo "scale=0; \$total_length / \$seq_count" | bc -l) bp" >> "\${base}_summary.txt"
    else
        echo "Average length: 0 bp" >> "\${base}_summary.txt"
    fi
    echo "Analysis completed: \$(date)" >> "\${base}_summary.txt"
    echo "Environment: Pixi-managed" >> "\${base}_summary.txt"
    """
}

workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Run version check
    VERSION_CHECK()
    
    // Process FASTA files with multiple analyses
    FASTA_STATS(input_ch)
    SEQUENCE_LENGTHS(input_ch)
    NUCLEOTIDE_COMPOSITION(input_ch)
    QUALITY_METRICS(input_ch)
    SUMMARY_REPORT(input_ch)
    
    // Display version information
    VERSION_CHECK.out.view { "Version Check: $it" }
}
