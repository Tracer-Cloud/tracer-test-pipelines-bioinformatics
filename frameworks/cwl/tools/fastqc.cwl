#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
baseCommand: fastqc

requirements:
  - class: InlineJavascriptRequirement

inputs:
  fastq_file:
    type: File
    inputBinding:
      position: 1
  outdir:
    type: string
    default: "."
    inputBinding:
      prefix: --outdir=
      separate: false

outputs:
  html_report:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file.basename.split('.')[0])_fastqc.html"
  zip_report:
    type: File
    outputBinding:
      glob: "$(inputs.fastq_file.basename.split('.')[0])_fastqc.zip"