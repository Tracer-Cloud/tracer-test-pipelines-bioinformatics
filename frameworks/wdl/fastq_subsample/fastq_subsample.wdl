version 1.0

workflow fastq_subsample_workflow {

  # Require first read, optional second read
  # Require a sample identifier (sample_id) to name the subsampled file
  input {
    File fastqgz_file_read_1
    File? fastqgz_file_read_2
    String sample_id
  }
  
  # See if a second read is present
  Boolean is_paired_end = defined(fastqgz_file_read_2)
  
  call sample_file {
    input: 
      fastqgz_file_read_1 = fastqgz_file_read_1,
      fastqgz_file_read_2 = fastqgz_file_read_2,
      sample_id = sample_id,
      is_paired_end = is_paired_end
  }
  
  output {
    File read1_subsample = sample_file.read1_subsample
    File? read2_subsample = sample_file.read2_subsample
  }

}

task sample_file {

  input {
    File fastqgz_file_read_1
    File? fastqgz_file_read_2
    String sample_id
    Int n = 10000
    Boolean is_paired_end
  }

  # Use seqtk to subsample files
  command <<<
  
    echo "Creating a subsampled file with ~{n} lines."
    seqtk sample -s11 ~{fastqgz_file_read_1} ~{n} > "~{sample_id}_1_subsample.fq"
    echo "Created ~{sample_id}_1_subsample.fq."
    
    if [ ~{is_paired_end} = true ] ; then
      echo "Subsampling paired-end reads for ~{sample_id}."
      seqtk sample -s11 ~{fastqgz_file_read_2} ~{n} > "~{sample_id}_2_subsample.fq"
      echo "Created ~{sample_id}_2_subsample.fq."
    fi
    
  >>>

  output {
    File read1_subsample="~{sample_id}_1_subsample.fq"
    File? read2_subsample="~{sample_id}_2_subsample.fq"
  }

  runtime {
    docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    memory:"8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible:  1
  }

}
