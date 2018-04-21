#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "STAR genomeGenerate" 
label: Spliced Transcripts Alignment to a Reference

requirements:
  DockerRequirement:
    dockerPull: welliton/star:v2.5.3a

baseCommand: [STAR, --runMode, genomeGenerate, --genomeDir, .]

inputs:

  genomeFastaFiles:
    type: File[]
    inputBinding:
      prefix: --genomeFastaFiles
    doc: Paths to the fasta files with genomic sequences for genome generation.
  sjdbGTFfile:
    type: File?
    inputBinding:
      prefix: --sjdbGTFfile
    doc: Path to the GTF file with annotations.

outputs:
  index_files:
    type: File[]
    outputBinding:
      glob: "*"
  