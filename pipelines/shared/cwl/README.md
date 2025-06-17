# How to run on Tracer Sandbox

## Requirements 
### 1. Install cwltool 
```bash
 # enable the universe repo
sudo add-apt-repository universe         

# refresh package lists
sudo apt update                       

# Install cwltool
sudo apt install cwltool

# Check if it works
cwltool --version
```


### 2. Creation of dummy files
- You need to create small input files that will be the input for the RNAseq job. 
- All dummy files can be recreated with 3 helper commands.  It keeps the repository clean by writing into integrations/cwl/.


What is important is to create (put them in the right folder)
- a dummy annotation file (test.gtf)
- a simple reference genome (genome.fa)
- tiny fake fastq file (test.fastq)

#### Test files to create 
##### Creation of dummy annotation file
```bash
cat > test.gtf <<EOF
chr1\ttest\texon\t1\t20\t.\t+\t.\tgene_id "gene1";
EOF
```

##### Make a toy reference genome
```bash
mkdir -p test-index
echo -e ">chr1\nACGTACGTACGTACGTACGT" > test-index/genome.fa
hisat2-build test-index/genome.fa test-index/genome
```

##### Creation of tiny fake fastq file
```bash
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
```

### 3. Run the CWL Command
### Run the command 
```bash
make run_cwl
```

### 4. Results you should see when the pipeline is working
![image](https://github.com/user-attachments/assets/bf391ef8-f056-4965-897c-d8d06200d4c3)


