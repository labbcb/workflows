#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "bowtie" 
label: "An ultrafast memory-efficient short read aligner."
doc: |
  Bowtie is an ultrafast, memory-efficient short read aligner. It aligns short 
  DNA sequences (reads) to the human genome at a rate of over 25 million 35-bp 
  reads per hour. Bowtie indexes the genome with a Burrows-Wheeler index to keep 
  its memory footprint small.

baseCommand: bowtie

requirements:
  InitialWorkDirRequirement:
    listing: $(inputs.index_files)
  DockerRequirement:
    dockerPull: welliton/bowtie:v1.2.0

inputs:

# Main arguments

  index_files:
    type: File[]
    doc: Index files (*.ebwt) created by bowtie-index.
  ebwt:
    type: string
    inputBinding:
      position: 2
    doc: The basename of the index to be searched.
  file:
    type: File
    inputBinding:
      position: 3
    doc: |
      A single file containing unpaired reads to be aligned.
  hit:
    type: string
    inputBinding:
      position: 4
    doc: |
      File to write alignments to.
      By default, alignments are written to the "standard out" filehandle.

# Alignment

  seedmms:
    type: int
    default: 2
    inputBinding:
      position: 1
      prefix: --seedmms
    doc: |
      Maximum number of mismatches permitted in the "seed".
      This may be 0, 1, 2 or 3 and the default is 2.

# Reporting
  all:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: --all
    doc: Report all valid alignments per read or pair.
  best:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: --best
    doc: |
      Make Bowtie guarantee that reported singleton alignments are "best"
      in terms of stratum.
  strata:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: --strata
    doc: |
      If many valid alignments exist and are reportable and they fall into more
      than one alignment "stratum", report only those alignments that fall into 
      the best stratum.

# SAM

  sam:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: --sam
    doc: Print alignments in SAM format.

# Performance

  threads:
    type: int
    default: 1
    inputBinding:
      position: 1
      prefix: --threads
    doc: "Launch parallel search threads (default: 1)."
    
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.hit)
