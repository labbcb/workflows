#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "htseq-count"
label: "Analysing high-throughput sequencing data with Python"

baseCommand: htseq-count

requirements:
   DockerRequirement:
    dockerPull: welliton/htseq:v0.9.1

inputs:

  input:
    type: File
    inputBinding:
      position: 2
    doc: Input file.
  gtf_file:
    type: File
    inputBinding:
      position: 3
    doc: GTF File.
  f:
    type: string
    default: sam
    inputBinding:
      prefix: --format
      position: 1
    doc: |
      Format of the input data. Possible values are sam (for text SAM files) 
      and bam (for binary BAM files). Default is sam.
  stranded:
    type: string
    default: yes
    inputBinding:
      prefix: --stranded
      position: 1
    doc: |
      For stranded=no, a read is considered overlapping with a feature 
      regardless of whether it is mapped to the same or the opposite strand as 
      the feature. For stranded=yes and single-end reads, the read has to be 
      mapped to the same strand as the feature. For paired-end reads, the first 
      read has to be on the same strand and the second read on the opposite 
      strand. For stranded=reverse, these rules are reversed.

stdout: "$(inputs.input.basename)_count.txt"

outputs:

  output:
    type: File
    outputBinding:
      glob: "$(inputs.input.basename)_count.txt"
