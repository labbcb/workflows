#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "STAR"
label: Spliced Transcripts Alignment to a Reference

requirements:
  InitialWorkDirRequirement:
    listing: $(inputs.index_files)
  DockerRequirement:
    dockerPull: welliton/star:2.5.3a

baseCommand: [STAR, --runMode, alignReads, --genomeDir, .]

inputs:
  index_files:
    type: File[]
  readFilesIn_R1:
    type: File
    inputBinding:
      prefix: --readFilesIn
      position: 1
  readFilesIn_R2:
    type: File
    inputBinding:
      position: 2
  outSAMtype:
    type: string[]?
    inputBinding:
      prefix: --outSAMtype
  outFileNamePrefix:
    type: string
    inputBinding:
      prefix: --outFileNamePrefix

outputs:
  output:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix)Aligned.*"
  log_final:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix)Log.final.out"
  log:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix)Log.out"
  log_progress:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix)Log.progress.out"
  sj:
    type: File
    outputBinding:
      glob: "$(inputs.outFileNamePrefix)SJ.out.tab"
