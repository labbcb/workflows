cwlVersion: v1.0
class: CommandLineTool
id: bwa-index
label: "BWA index"
doc: Index sequences in the FASTA format

baseCommand: [bwa, index]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.file) ]

hints:
  DockerRequirement:
    dockerPull: welliton/bwa:0.7.15

inputs:
  file:
    type: File
    inputBinding:
      position: 1
  a:
    type: string?
    doc: BWT construction algorithm bwtsw, is or rb2 [auto]
    inputBinding:
      prefix: -a
  p:
    type: string?
    doc: prefix of the index [same as fasta name]
    inputBinding:
      prefix: -p
  b:
    type: int?
    doc: block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
    inputBinding:
      prefix: -b
  sixtyfour:
    type: boolean
    default: false
    doc: index files named as <in.fasta>.64.* instead of <in.fasta>.*
    inputBinding:
      prefix: "-6"

outputs:
  output:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    outputBinding:
      glob: $(inputs.file.basename)
