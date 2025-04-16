# How to run on Tracer Sandbox
```bash
cwltool rnaseq-workflow.cwl  rnaseq-job-001.yml 
```

todo's to run a minimal viable prototype.
- You need to have a small input file such as SRR3907288_1.fastq this will be an input for the RNAseq job. 
- You need a tiny test index. 

# Creation of dummy annotation file
cat > test.gtf <<EOF
chr1\ttest\texon\t1\t20\t.\t+\t.\tgene_id "gene1";
EOF

# Make a toy reference genome
mkdir -p test-index
echo -e ">chr1\nACGTACGTACGTACGTACGT" > test-index/genome.fa
hisat2-build test-index/genome.fa test-index/genome


# Creation of tiny fake fastq file
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
