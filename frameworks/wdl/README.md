# Tracer WDL exmaple
- Istall Java: 

```
bashsudo apt-get update && sudo apt-get install default-jre -y  
```

- Install Cromwell: 

```bash
wget https://github.com/broadinstitute/cromwell/releases/download/84/cromwell-84.jar -O cromwell.jar
```
```bash
# Verify installation
java -jar cromwell.jar --version
```


fastq_subsample.wdl)

# Other notes

WDL workflows for use on AnVIL and other platforms.

- `fastq_subsample`: subsets, or samples, a fastq.gz file. It currently takes a number of random reads to create a smaller file. This is intended for the purposes of educational activities, creating files for testing, etc. [Dockstore link](https://dockstore.org/workflows/github.com/fhdsl/AnVIL_WDLs/fastq_subsample).



# Create a simple fastqfile 