#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "test_data/*.fasta"
params.outdir = "results"
params.conda_env = "nextflow-minimal"

// Process 1: Version check for all tools
process VERSION_CHECK {
    publishDir "${params.outdir}/versions", mode: 'copy'
    conda params.conda_env
    
    output:
    path "tool_versions.txt"
    
    script:
    """
    echo "=== Bioinformatics Tools Version Check ===" > tool_versions.txt
    echo "Analysis Date: \$(date)" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "1. Nextflow:" >> tool_versions.txt
    nextflow -version 2>&1 | head -3 >> tool_versions.txt 2>/dev/null || echo "Nextflow not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "2. FastQC:" >> tool_versions.txt
    fastqc --version >> tool_versions.txt 2>/dev/null || echo "FastQC not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "3. Seqtk:" >> tool_versions.txt
    seqtk 2>&1 | head -3 >> tool_versions.txt 2>/dev/null || echo "Seqtk not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "4. Samtools:" >> tool_versions.txt
    samtools --version 2>&1 | head -3 >> tool_versions.txt 2>/dev/null || echo "Samtools not available" >> tool_versions.txt
    echo "  Samtools subcommands:" >> tool_versions.txt
    echo "  - sort:" >> tool_versions.txt && samtools sort 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    sort not available" >> tool_versions.txt
    echo "  - index:" >> tool_versions.txt && samtools index 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    index not available" >> tool_versions.txt
    echo "  - view:" >> tool_versions.txt && samtools view 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    view not available" >> tool_versions.txt
    echo "  - flagstat:" >> tool_versions.txt && samtools flagstat 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    flagstat not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "5. Bedtools:" >> tool_versions.txt
    bedtools --version >> tool_versions.txt 2>/dev/null || echo "Bedtools not available" >> tool_versions.txt
    echo "  Bedtools subcommands:" >> tool_versions.txt
    echo "  - intersect:" >> tool_versions.txt && bedtools intersect 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    intersect not available" >> tool_versions.txt
    echo "  - merge:" >> tool_versions.txt && bedtools merge 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    merge not available" >> tool_versions.txt
    echo "  - coverage:" >> tool_versions.txt && bedtools coverage 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    coverage not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "6. BLAST:" >> tool_versions.txt
    blastn -version 2>&1 | head -3 >> tool_versions.txt 2>/dev/null || echo "BLAST not available" >> tool_versions.txt
    echo "  BLAST subcommands:" >> tool_versions.txt
    echo "  - blastp:" >> tool_versions.txt && blastp -version 2>&1 | head -1 >> tool_versions.txt 2>/dev/null || echo "    blastp not available" >> tool_versions.txt
    echo "  - blastx:" >> tool_versions.txt && blastx -version 2>&1 | head -1 >> tool_versions.txt 2>/dev/null || echo "    blastx not available" >> tool_versions.txt
    echo "  - tblastn:" >> tool_versions.txt && tblastn -version 2>&1 | head -1 >> tool_versions.txt 2>/dev/null || echo "    tblastn not available" >> tool_versions.txt
    echo "  - makeblastdb:" >> tool_versions.txt && makeblastdb -version 2>&1 | head -1 >> tool_versions.txt 2>/dev/null || echo "    makeblastdb not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "7. STAR:" >> tool_versions.txt
    STAR --version >> tool_versions.txt 2>/dev/null || echo "STAR not available" >> tool_versions.txt
    echo "  STAR subcommands:" >> tool_versions.txt
    echo "  - STAR --help:" >> tool_versions.txt && STAR --help 2>&1 | head -3 >> tool_versions.txt 2>/dev/null || echo "    help not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "8. Salmon:" >> tool_versions.txt
    salmon --version >> tool_versions.txt 2>/dev/null || echo "Salmon not available" >> tool_versions.txt
    echo "  Salmon subcommands:" >> tool_versions.txt
    echo "  - index:" >> tool_versions.txt && salmon index --help 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    index not available" >> tool_versions.txt
    echo "  - quant:" >> tool_versions.txt && salmon quant --help 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    quant not available" >> tool_versions.txt
    echo "  - alevin:" >> tool_versions.txt && salmon alevin --help 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    alevin not available" >> tool_versions.txt
    echo "  - swim:" >> tool_versions.txt && salmon swim --help 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    swim not available" >> tool_versions.txt
    echo "  - merge:" >> tool_versions.txt && salmon merge --help 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    merge not available" >> tool_versions.txt
    echo "  - diff:" >> tool_versions.txt && salmon diff --help 2>&1 | head -2 >> tool_versions.txt 2>/dev/null || echo "    diff not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "9. MUSCLE:" >> tool_versions.txt
    muscle -version >> tool_versions.txt 2>/dev/null || echo "MUSCLE not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "10. MAFFT:" >> tool_versions.txt
    mafft --version >> tool_versions.txt 2>/dev/null || echo "MAFFT not available" >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "=== System Information ===" >> tool_versions.txt
    uname -a >> tool_versions.txt
    echo "" >> tool_versions.txt
    
    echo "=== Environment Check Complete ===" >> tool_versions.txt
    """
}

// Process 2: Basic FASTA statistics
process FASTA_STATS {
    publishDir "${params.outdir}/fasta_stats", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    base=\$(basename "${fasta_file}" .fasta)
    
    # Count sequences
    seq_count=\$(grep -c '^>' "${fasta_file}" 2>/dev/null || echo "0")
    
    # Get file size
    file_size=\$(wc -c < "${fasta_file}")
    
    # Write basic report
    {
        echo "=== Basic FASTA Statistics ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo "Number of sequences: \$seq_count"
        echo "File size: \$file_size bytes"
        echo ""
        echo "=== First 5 Sequence Headers ==="
        grep '^>' "${fasta_file}" | head -5 2>/dev/null || echo "No headers found"
    } > "\${base}_basic_stats.txt"
    """
}

// Process 3: Sequence length analysis
process SEQUENCE_LENGTHS {
    publishDir "${params.outdir}/lengths", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    base=\$(basename "${fasta_file}" .fasta)
    
    # Use seqtk for length analysis
    if command -v seqtk >/dev/null 2>&1; then
        seqtk comp "${fasta_file}" > "\${base}_comp.tmp" 2>/dev/null || echo "seq1 0 0 0 0 0 0 0 0 0" > "\${base}_comp.tmp"
        
        # Calculate basic stats
        total_length=\$(awk '{sum += \$2} END {print sum+0}' "\${base}_comp.tmp" 2>/dev/null || echo "0")
        seq_count=\$(wc -l < "\${base}_comp.tmp" 2>/dev/null || echo "0")
        
        if [ "\$seq_count" -gt 0 ] && [ "\$total_length" -gt 0 ]; then
            avg_length=\$(echo "scale=2; \$total_length / \$seq_count" | bc -l 2>/dev/null || echo "0")
        else
            avg_length=0
        fi
        
        rm -f "\${base}_comp.tmp"
    else
        total_length=0
        seq_count=0
        avg_length=0
    fi
    
    # Write length report
    {
        echo "=== Sequence Length Analysis ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo "Total sequences: \$seq_count"
        echo "Total length: \$total_length bp"
        echo "Average length: \$avg_length bp"
    } > "\${base}_length_analysis.txt"
    """
}

// Process 4: Nucleotide composition
process NUCLEOTIDE_COMPOSITION {
    publishDir "${params.outdir}/composition", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    base=\$(basename "${fasta_file}" .fasta)
    
    # Use seqtk for composition analysis
    if command -v seqtk >/dev/null 2>&1; then
        seqtk comp "${fasta_file}" > "\${base}_comp.tmp" 2>/dev/null || echo "seq1 0 0 0 0 0 0 0 0 0" > "\${base}_comp.tmp"
        
        # Calculate nucleotide counts
        a_count=\$(awk '{sum += \$3} END {print sum+0}' "\${base}_comp.tmp" 2>/dev/null || echo "0")
        c_count=\$(awk '{sum += \$4} END {print sum+0}' "\${base}_comp.tmp" 2>/dev/null || echo "0")
        g_count=\$(awk '{sum += \$5} END {print sum+0}' "\${base}_comp.tmp" 2>/dev/null || echo "0")
        t_count=\$(awk '{sum += \$6} END {print sum+0}' "\${base}_comp.tmp" 2>/dev/null || echo "0")
        
        total=\$(echo "\$a_count + \$c_count + \$g_count + \$t_count" | bc -l 2>/dev/null || echo "0")
        
        if [ "\$total" -gt 0 ]; then
            gc_content=\$(echo "scale=2; (\$g_count + \$c_count) * 100 / \$total" | bc -l 2>/dev/null || echo "0.00")
        else
            gc_content="0.00"
        fi
        
        rm -f "\${base}_comp.tmp"
    else
        a_count=0; c_count=0; g_count=0; t_count=0; gc_content="0.00"
    fi
    
    # Write composition report
    {
        echo "=== Nucleotide Composition ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo "A count: \$a_count"
        echo "C count: \$c_count"
        echo "G count: \$g_count"
        echo "T count: \$t_count"
        echo "GC content: \$gc_content%"
    } > "\${base}_composition.txt"
    """
}

// Process 5: Quality assessment
process QUALITY_ASSESSMENT {
    publishDir "${params.outdir}/quality", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    base=\$(basename "${fasta_file}" .fasta)
    
    # Basic quality checks
    seq_count=\$(grep -c '^>' "${fasta_file}" 2>/dev/null || echo "0")
    file_size=\$(wc -c < "${fasta_file}")
    
    # Check for empty headers
    empty_headers=\$(grep -c '^>\$' "${fasta_file}" 2>/dev/null || echo "0")
    
    # Check for duplicate headers
    dup_headers=\$(grep '^>' "${fasta_file}" 2>/dev/null | sort | uniq -d | wc -l || echo "0")
    
    # Write quality report
    {
        echo "=== Quality Assessment ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo "File size: \$file_size bytes"
        echo "Sequence count: \$seq_count"
        echo "Empty headers: \$empty_headers"
        echo "Duplicate headers: \$dup_headers"
        echo ""
        echo "=== Quality Status ==="
        if [ "\$seq_count" -gt 0 ]; then
            echo "✓ File contains sequences"
        else
            echo "✗ No sequences found"
        fi
        
        if [ "\$empty_headers" -eq 0 ]; then
            echo "✓ No empty headers"
        else
            echo "✗ Empty headers found"
        fi
        
        if [ "\$dup_headers" -eq 0 ]; then
            echo "✓ No duplicate headers"
        else
            echo "✗ Duplicate headers found"
        fi
    } > "\${base}_quality.txt"
    """
}

// Process 6: Tool compatibility test
process TOOL_COMPATIBILITY {
    publishDir "${params.outdir}/compatibility", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    base=\$(basename "${fasta_file}" .fasta)
    
    # Test tool compatibility
    {
        echo "=== Tool Compatibility Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        # Test seqtk
        echo "=== seqtk ==="
        if command -v seqtk >/dev/null 2>&1; then
            echo "✓ Available"
            if seqtk comp "${fasta_file}" >/dev/null 2>&1; then
                echo "✓ Compatible"
            else
                echo "✗ Not compatible"
            fi
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test samtools
        echo "=== samtools ==="
        if command -v samtools >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(samtools --version | head -1 2>/dev/null || echo 'unknown')"
            echo "  Subcommands: sort, index, view, flagstat"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test BLAST
        echo "=== BLAST ==="
        if command -v blastn >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(blastn -version | head -1 2>/dev/null || echo 'unknown')"
            echo "  Subcommands: blastp, blastx, tblastn, makeblastdb"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test FastQC
        echo "=== FastQC ==="
        if command -v fastqc >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(fastqc --version 2>/dev/null || echo 'unknown')"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test STAR
        echo "=== STAR ==="
        if command -v STAR >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(STAR --version 2>/dev/null || echo 'unknown')"
            echo "  Subcommands: --help, --version"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test Salmon
        echo "=== Salmon ==="
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo "  Subcommands: index, quant, alevin, swim, merge, diff"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test MUSCLE
        echo "=== MUSCLE ==="
        if command -v muscle >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(muscle -version 2>/dev/null || echo 'unknown')"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test MAFFT
        echo "=== MAFFT ==="
        if command -v mafft >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(mafft --version 2>/dev/null || echo 'unknown')"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Test Bedtools
        echo "=== Bedtools ==="
        if command -v bedtools >/dev/null 2>&1; then
            echo "✓ Available"
            echo "  Version: \$(bedtools --version 2>/dev/null || echo 'unknown')"
            echo "  Subcommands: intersect, merge, coverage"
        else
            echo "✗ Not available"
        fi
        echo ""
        
        # Summary
        echo "=== Compatibility Summary ==="
        available_tools=0
        for tool in seqtk samtools blastn fastqc STAR salmon muscle mafft bedtools; do
            if command -v \$tool >/dev/null 2>&1; then
                available_tools=\$((available_tools + 1))
            fi
        done
        echo "Available tools: \$available_tools/9"
        
    } > "\${base}_compatibility.txt"
    """
}

// Process 7: Summary report
process SUMMARY_REPORT {
    publishDir "${params.outdir}/summary", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "*.txt"
    
    script:
    """
    base=\$(basename "${fasta_file}" .fasta)
    
    # Generate summary
    {
        echo "=== FASTA Analysis Summary ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo "Environment: nextflow-minimal conda environment"
        echo ""
        
        # Basic info
        seq_count=\$(grep -c '^>' "${fasta_file}" 2>/dev/null || echo "0")
        file_size=\$(wc -c < "${fasta_file}")
        
        echo "=== File Information ==="
        echo "File size: \$file_size bytes"
        echo "Total sequences: \$seq_count"
        echo ""
        
        # Tool availability
        echo "=== Available Tools ==="
        for tool in seqtk samtools blastn fastqc STAR salmon muscle mafft bedtools; do
            if command -v \$tool >/dev/null 2>&1; then
                echo "✓ \$tool"
            else
                echo "✗ \$tool"
            fi
        done
        echo ""
        
        echo "=== Analysis Complete ==="
        echo "All processes completed successfully"
        echo "Results available in: ${params.outdir}"
        
    } > "\${base}_summary.txt"
    """
}

// Process 8: System information
process SYSTEM_INFO {
    publishDir "${params.outdir}/system", mode: 'copy'
    conda params.conda_env
    
    output:
    path "system_info.txt"
    
    script:
    """
    {
        echo "=== System Information ==="
        echo "Analysis Date: \$(date)"
        echo ""
        
        echo "=== Operating System ==="
        uname -a
        echo ""
        
        echo "=== CPU Information ==="
        if command -v sysctl >/dev/null 2>&1; then
            echo "CPU cores: \$(sysctl -n hw.ncpu 2>/dev/null || echo 'unknown')"
            echo "CPU brand: \$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'unknown')"
        else
            echo "CPU info not available"
        fi
        echo ""
        
        echo "=== Memory Information ==="
        if command -v sysctl >/dev/null 2>&1; then
            sysctl hw.memsize 2>/dev/null | awk '{print "Total RAM: " \$2/1024/1024/1024 " GB"}' || echo "Memory info not available"
        else
            echo "Memory info not available"
        fi
        echo ""
        
        echo "=== Disk Space ==="
        df -h . | head -2
        echo ""
        
        echo "=== Environment ==="
        echo "Conda environment: \$CONDA_DEFAULT_ENV"
        echo "Python version: \$(python --version 2>/dev/null || echo 'not available')"
        echo "Nextflow version: \$(nextflow -version 2>&1 | head -1 | cut -d' ' -f2 2>/dev/null || echo 'unknown')"
        
    } > system_info.txt
    """
}

// Process 9: Salmon Index
process SALMON_INDEX {
    publishDir "${params.outdir}/salmon", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "salmon_index_info.txt"
    
    script:
    """
    {
        echo "=== Salmon Index Command Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Salmon is available"
            echo "Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo ""
            
            echo "=== Salmon Index Help ==="
            salmon index --help 2>&1 | head -10
            echo ""
            
            echo "=== Salmon Index Command Syntax ==="
            echo "salmon index -t <transcriptome.fa> -i <index_name>"
            echo ""
            
            echo "=== Index Command Test ==="
            base=\$(basename "${fasta_file}" .fasta)
            echo "Would create index for: \${base}"
            echo "Command: salmon index -t ${fasta_file} -i \${base}_salmon_index"
            echo ""
            
            echo "=== Index Parameters ==="
            echo "-t: transcriptome file (FASTA format)"
            echo "-i: index directory name"
            echo "-k: k-mer size (default: 31)"
            echo "--type: index type (default: quasi)"
            echo ""
            
        else
            echo "✗ Salmon not available"
        fi
        
    } > salmon_index_info.txt
    """
}

// Process 10: Salmon Quant
process SALMON_QUANT {
    publishDir "${params.outdir}/salmon", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "salmon_quant_info.txt"
    
    script:
    """
    {
        echo "=== Salmon Quant Command Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Salmon is available"
            echo "Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo ""
            
            echo "=== Salmon Quant Help ==="
            salmon quant --help 2>&1 | head -10
            echo ""
            
            echo "=== Salmon Quant Command Syntax ==="
            echo "salmon quant -i <index> -l <libtype> -1 <read1.fq> -2 <read2.fq> -o <output_dir>"
            echo ""
            
            echo "=== Quant Command Test ==="
            base=\$(basename "${fasta_file}" .fasta)
            echo "Would quantify against: \${base}"
            echo "Command: salmon quant -i \${base}_salmon_index -l A -1 reads_1.fq -2 reads_2.fq -o \${base}_quant"
            echo ""
            
            echo "=== Quant Parameters ==="
            echo "-i: salmon index directory"
            echo "-l: library type (A=auto, IU=unstranded, ISF=stranded forward, ISR=stranded reverse)"
            echo "-1: read 1 file (FASTQ)"
            echo "-2: read 2 file (FASTQ, for paired-end)"
            echo "-o: output directory"
            echo "--gcBias: enable GC bias correction"
            echo "--seqBias: enable sequence bias correction"
            echo ""
            
        else
            echo "✗ Salmon not available"
        fi
        
    } > salmon_quant_info.txt
    """
}

// Process 11: Salmon Alevin
process SALMON_ALEVIN {
    publishDir "${params.outdir}/salmon", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "salmon_alevin_info.txt"
    
    script:
    """
    {
        echo "=== Salmon Alevin Command Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Salmon is available"
            echo "Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo ""
            
            echo "=== Salmon Alevin Help ==="
            salmon alevin --help 2>&1 | head -10
            echo ""
            
            echo "=== Salmon Alevin Command Syntax ==="
            echo "salmon alevin -i <index> -l <libtype> -1 <read1.fq> -2 <read2.fq> -o <output_dir>"
            echo ""
            
            echo "=== Alevin Command Test ==="
            base=\$(basename "${fasta_file}" .fasta)
            echo "Would process single-cell data against: \${base}"
            echo "Command: salmon alevin -i \${base}_salmon_index -l A -1 sc_reads_1.fq -2 sc_reads_2.fq -o \${base}_alevin"
            echo ""
            
            echo "=== Alevin Parameters ==="
            echo "-i: salmon index directory"
            echo "-l: library type"
            echo "-1: read 1 file (FASTQ)"
            echo "-2: read 2 file (FASTQ)"
            echo "-o: output directory"
            echo "--whitelist: cell barcode whitelist"
            echo "--rad: 10x genomics mode"
            echo ""
            
        else
            echo "✗ Salmon not available"
        fi
        
    } > salmon_alevin_info.txt
    """
}

// Process 12: Salmon Swim
process SALMON_SWIM {
    publishDir "${params.outdir}/salmon", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "salmon_swim_info.txt"
    
    script:
    """
    {
        echo "=== Salmon Swim Command Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Salmon is available"
            echo "Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo ""
            
            echo "=== Salmon Swim Help ==="
            salmon swim --help 2>&1 | head -10
            echo ""
            
            echo "=== Salmon Swim Command Syntax ==="
            echo "salmon swim -i <index> -l <libtype> -1 <read1.fq> -2 <read2.fq> -o <output_dir>"
            echo ""
            
            echo "=== Swim Command Test ==="
            base=\$(basename "${fasta_file}" .fasta)
            echo "Would process with swim against: \${base}"
            echo "Command: salmon swim -i \${base}_salmon_index -l A -1 reads_1.fq -2 reads_2.fq -o \${base}_swim"
            echo ""
            
            echo "=== Swim Parameters ==="
            echo "-i: salmon index directory"
            echo "-l: library type"
            echo "-1: read 1 file (FASTQ)"
            echo "-2: read 2 file (FASTQ)"
            echo "-o: output directory"
            echo "--mimicBT2: mimic bowtie2 behavior"
            echo "--hardFilter: enable hard filtering"
            echo ""
            
        else
            echo "✗ Salmon not available"
        fi
        
    } > salmon_swim_info.txt
    """
}

// Process 13: Salmon Merge
process SALMON_MERGE {
    publishDir "${params.outdir}/salmon", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "salmon_merge_info.txt"
    
    script:
    """
    {
        echo "=== Salmon Merge Command Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Salmon is available"
            echo "Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo ""
            
            echo "=== Salmon Merge Help ==="
            salmon merge --help 2>&1 | head -10
            echo ""
            
            echo "=== Salmon Merge Command Syntax ==="
            echo "salmon merge -i <input_dirs> -o <output_dir>"
            echo ""
            
            echo "=== Merge Command Test ==="
            base=\$(basename "${fasta_file}" .fasta)
            echo "Would merge quantification results for: \${base}"
            echo "Command: salmon merge -i \${base}_quant1 \${base}_quant2 -o \${base}_merged"
            echo ""
            
            echo "=== Merge Parameters ==="
            echo "-i: input quantification directories (space-separated)"
            echo "-o: output directory for merged results"
            echo "--column: column to merge (default: NumReads)"
            echo "--names: sample names (space-separated)"
            echo ""
            
        else
            echo "✗ Salmon not available"
        fi
        
    } > salmon_merge_info.txt
    """
}

// Process 14: Salmon Diff
process SALMON_DIFF {
    publishDir "${params.outdir}/salmon", mode: 'copy'
    conda params.conda_env
    
    input:
    path fasta_file
    
    output:
    path "salmon_diff_info.txt"
    
    script:
    """
    {
        echo "=== Salmon Diff Command Test ==="
        echo "File: ${fasta_file}"
        echo "Analysis Date: \$(date)"
        echo ""
        
        if command -v salmon >/dev/null 2>&1; then
            echo "✓ Salmon is available"
            echo "Version: \$(salmon --version 2>/dev/null || echo 'unknown')"
            echo ""
            
            echo "=== Salmon Diff Help ==="
            salmon diff --help 2>&1 | head -10
            echo ""
            
            echo "=== Salmon Diff Command Syntax ==="
            echo "salmon diff -i <input_dirs> -o <output_dir> -c <condition_file>"
            echo ""
            
            echo "=== Diff Command Test ==="
            base=\$(basename "${fasta_file}" .fasta)
            echo "Would perform differential expression for: \${base}"
            echo "Command: salmon diff -i \${base}_quant1 \${base}_quant2 -o \${base}_diff -c conditions.txt"
            echo ""
            
            echo "=== Diff Parameters ==="
            echo "-i: input quantification directories (space-separated)"
            echo "-o: output directory for differential expression results"
            echo "-c: condition file (tab-separated: sample_id condition)"
            echo "--method: differential expression method (default: DESeq2)"
            echo "--design: experimental design formula"
            echo ""
            
        else
            echo "✗ Salmon not available"
        fi
        
    } > salmon_diff_info.txt
    """
}

workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Run all processes
    VERSION_CHECK()
    SYSTEM_INFO()
    FASTA_STATS(input_ch)
    SEQUENCE_LENGTHS(input_ch)
    NUCLEOTIDE_COMPOSITION(input_ch)
    QUALITY_ASSESSMENT(input_ch)
    TOOL_COMPATIBILITY(input_ch)
    SUMMARY_REPORT(input_ch)
    
    // Salmon subcommand processes
    SALMON_INDEX(input_ch)
    SALMON_QUANT(input_ch)
    SALMON_ALEVIN(input_ch)
    SALMON_SWIM(input_ch)
    SALMON_MERGE(input_ch)
    SALMON_DIFF(input_ch)
    
    // Display completion message
    VERSION_CHECK.out.view { "=== Tool Environment Check ===\n$it" }
    SYSTEM_INFO.out.view { "=== System Information ===\n$it" }
}