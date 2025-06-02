#!/usr/bin/env bash
set -euo pipefail

# Create test.gtf
cat > test.gtf <<EOF
chr1	test	exon	1	20	.	+	.	gene_id "gene1";
EOF

echo "✅ Created test.gtf"

# Create reference genome fasta and index it
mkdir -p test-index
cat > test-index/genome.fa <<EOF
>chr1
ACGTACGTACGTACGTACGT
EOF

echo "✅ Created test-index/genome.fa"

# Build index using hisat2-build
if command -v hisat2-build &> /dev/null; then
    hisat2-build test-index/genome.fa test-index/genome
    echo "✅ Built hisat2 index"
else
    echo "❌ hisat2-build not found. Please install HISAT2 and retry."
    exit 1
fi

# Create small FASTQ file
cat > test.fastq <<EOF
@read1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read3
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
EOF

echo "✅ Created test.fastq"
