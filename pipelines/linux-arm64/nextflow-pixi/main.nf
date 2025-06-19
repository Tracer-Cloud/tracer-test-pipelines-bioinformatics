#!/usr/bin/env nextflow

// Simple pipeline that processes FASTA files
// Uses FastQC and STAR for version checking and basic analysis
// Optimized for Linux ARM64 with Pixi environment management

nextflow.enable.dsl = 2

// Parameters
params.input = "test_data/*.fasta"
params.outdir = "results"
params.help = false

// Help message
def helpMessage() {
    log.info"""
    ===================================
    Nextflow Pixi Pipeline (Linux ARM64)
    ===================================
    
    Usage:
    nextflow run main.nf [options]
    
    Options:
    --input         Input FASTA files (default: "test_data/*.fasta")
    --outdir        Output directory (default: "results")
    --help          Show this help message
    
    Example:
    nextflow run main.nf --input "data/*.fasta" --outdir "my_results"
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Process 1: Check FastQC and STAR versions
process CHECK_VERSIONS {
    label 'version_check'
    publishDir "${params.outdir}/versions", mode: 'copy'

    output:
    path "tool_versions.txt"

    script:
    """
    echo "=== Tool Versions in Pixi Environment (Linux ARM64) ===" > tool_versions.txt
    echo "Generated on: \$(date)" >> tool_versions.txt
    echo "Environment: Linux ARM64 with Pixi" >> tool_versions.txt
    echo "Architecture: \$(uname -m)" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "=== FastQC Version ===" >> tool_versions.txt
    fastqc --version >> tool_versions.txt 2>&1
    echo "" >> tool_versions.txt
    
    echo "=== STAR Version ===" >> tool_versions.txt
    STAR --version >> tool_versions.txt 2>&1
    echo "" >> tool_versions.txt
    
    echo "=== System Information ===" >> tool_versions.txt
    echo "OS: \$(uname -a)" >> tool_versions.txt
    echo "CPU cores: \$(nproc)" >> tool_versions.txt
    echo "Memory: \$(free -h | grep '^Mem:' | awk '{print \$2}')" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "=== Pixi Environment ===" >> tool_versions.txt
    echo "Pixi managed environment active" >> tool_versions.txt
    echo "All tools available in PATH" >> tool_versions.txt
    """
}

// Process 2: Run FastQC on FASTA files
process FASTQC_ANALYSIS {
    tag "$fasta_file.baseName"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_fastqc.html", optional: true
    path "${fasta_file.baseName}_fastqc.zip", optional: true
    path "${fasta_file.baseName}_fastqc_summary.txt"

    script:
    """
    echo "=== FastQC Analysis for ${fasta_file} ===" > ${fasta_file.baseName}_fastqc_summary.txt
    echo "File: ${fasta_file}" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "Generated on: \$(date)" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "Environment: Linux ARM64 with Pixi" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "" >> ${fasta_file.baseName}_fastqc_summary.txt
    
    # Run FastQC (note: FastQC expects FASTQ files, but we'll run it anyway for demonstration)
    echo "Running FastQC analysis..." >> ${fasta_file.baseName}_fastqc_summary.txt
    
    # FastQC might not work perfectly with FASTA files, so we'll capture the attempt
    if fastqc ${fasta_file} --outdir . 2>&1 | tee -a ${fasta_file.baseName}_fastqc_summary.txt; then
        echo "FastQC completed successfully" >> ${fasta_file.baseName}_fastqc_summary.txt
    else
        echo "FastQC completed with warnings (expected for FASTA input)" >> ${fasta_file.baseName}_fastqc_summary.txt
    fi
    
    # Basic file statistics as fallback
    echo "" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "=== Basic File Statistics ===" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "File size: \$(du -h ${fasta_file} | cut -f1)" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "Number of sequences: \$(grep -c '^>' ${fasta_file} || echo '0')" >> ${fasta_file.baseName}_fastqc_summary.txt
    echo "Number of lines: \$(wc -l < ${fasta_file})" >> ${fasta_file.baseName}_fastqc_summary.txt
    """
}

// Process 3: STAR version and basic information
process STAR_INFO {
    publishDir "${params.outdir}/star", mode: 'copy'

    output:
    path "star_info.txt"

    script:
    """
    echo "=== STAR Information (Linux ARM64) ===" > star_info.txt
    echo "Generated on: \$(date)" >> star_info.txt
    echo "Environment: Linux ARM64 with Pixi" >> star_info.txt
    echo "" >> star_info.txt
    
    echo "=== STAR Version Details ===" >> star_info.txt
    STAR --version >> star_info.txt 2>&1
    echo "" >> star_info.txt
    
    echo "=== STAR Help (first 20 lines) ===" >> star_info.txt
    STAR --help 2>&1 | head -20 >> star_info.txt
    echo "... (truncated)" >> star_info.txt
    echo "" >> star_info.txt
    
    echo "=== STAR Parameters ===" >> star_info.txt
    echo "STAR is available and ready for genome alignment tasks" >> star_info.txt
    echo "This is a minimal example - full STAR usage requires genome indices" >> star_info.txt
    """
}

// Process 4: Basic FASTA statistics
process FASTA_STATS {
    tag "$fasta_file.baseName"
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
    echo "Environment: Linux ARM64 with Pixi" >> ${fasta_file.baseName}_stats.txt
    echo "" >> ${fasta_file.baseName}_stats.txt

    # Count sequences
    seq_count=\$(grep -c '^>' ${fasta_file} || echo "0")
    echo "Number of sequences: \$seq_count" >> ${fasta_file.baseName}_stats.txt

    # Calculate total length
    total_length=\$(grep -v '^>' ${fasta_file} | tr -d '\\n' | wc -c || echo "0")
    echo "Total sequence length: \$total_length bp" >> ${fasta_file.baseName}_stats.txt

    # Calculate average length
    if [ \$seq_count -gt 0 ]; then
        avg_length=\$(echo "scale=2; \$total_length / \$seq_count" | bc -l 2>/dev/null || echo "0")
        echo "Average sequence length: \$avg_length bp" >> ${fasta_file.baseName}_stats.txt
    else
        echo "Average sequence length: 0 bp" >> ${fasta_file.baseName}_stats.txt
    fi

    # Show sequence IDs (first 10)
    echo "" >> ${fasta_file.baseName}_stats.txt
    echo "Sequence IDs (first 10):" >> ${fasta_file.baseName}_stats.txt
    grep '^>' ${fasta_file} | head -10 >> ${fasta_file.baseName}_stats.txt
    
    if [ \$(grep -c '^>' ${fasta_file}) -gt 10 ]; then
        echo "... and \$((\$(grep -c '^>' ${fasta_file}) - 10)) more sequences" >> ${fasta_file.baseName}_stats.txt
    fi
    
    echo "" >> ${fasta_file.baseName}_stats.txt
    echo "Analysis completed with Pixi-managed tools on Linux ARM64" >> ${fasta_file.baseName}_stats.txt
    """
}

// Main workflow
workflow {
    // Show pipeline info
    log.info """
    ===================================
    Nextflow Pixi Pipeline (Linux ARM64)
    ===================================
    Input files: ${params.input}
    Output directory: ${params.outdir}
    Architecture: Linux ARM64
    Package manager: Pixi
    ===================================
    """
    
    // Create channel from input files
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Log the files being processed
    input_ch.view { "Processing file: $it" }
    
    // Run version checks
    CHECK_VERSIONS()
    STAR_INFO()
    
    // Run analysis processes
    FASTQC_ANALYSIS(input_ch)
    FASTA_STATS(input_ch)
    
    // Collect outputs
    fastqc_results = FASTQC_ANALYSIS.out[2].collect()
    stats_results = FASTA_STATS.out.collect()
    
    // Log completion
    fastqc_results.view { "Generated FastQC summaries: $it" }
    stats_results.view { "Generated statistics files: $it" }
}
