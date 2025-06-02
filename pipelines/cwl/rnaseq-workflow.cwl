cwlVersion: v1.2
class: Workflow

inputs:
  reads: File
  annotation: File
  hisat2_index: string
  feature_type:
    type: string
    default: "exon"
  output_file:
    type: string
    default: "counts.txt"

steps:
  qc:
    run: tools/fastqc.cwl
    in:
      fastq_file: reads
    out: [html_report]

  align:
    run: tools/hisat2.cwl
    in:
      reads: reads
      index_prefix: hisat2_index
    out: [sam_output]

  count:
    run: tools/featurecounts.cwl
    in:
      annotation: annotation
      sam_file: align/sam_output
      output_file: output_file
      feature_type: feature_type
    out: [count_output, summary_output]

outputs:
  report:
    type: ["null", File]
    outputSource: qc/html_report
    
  aligned_sam:
    type: File
    outputSource: align/sam_output

  counts:
    type: File
    outputSource: count/count_output
    
  count_summary:
    type: File
    outputSource: count/summary_output