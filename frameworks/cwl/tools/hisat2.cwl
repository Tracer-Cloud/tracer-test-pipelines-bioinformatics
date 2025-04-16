#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
baseCommand: hisat2

requirements:
  - class: InlineJavascriptRequirement

inputs:
  reads:
    type: File
    inputBinding:
      prefix: -U
      position: 2
  index_prefix:
    type: string  # This matches what the workflow is providing
    inputBinding:
      prefix: -x
      position: 1
  threads:
    type: int
    default: 6
    inputBinding:
      prefix: -p
      position: 3

stdout: output.sam

outputs:
  sam_output:
    type: stdout