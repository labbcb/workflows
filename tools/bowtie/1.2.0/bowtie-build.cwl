cwlVersion: v1.0
class: CommandLineTool

baseCommand: [bowtie-build]

requirements:
  DockerRequirement:
      dockerPull: welliton/bowtie:v1.2.0

inputs:
  genome_files:
    type: File[]
    doc: |
      A comma-separated list of FASTA files containing the reference sequences
      to be aligned to.
    inputBinding:
      position: 1
      itemSeparator: ","
  ebwt_base:
    type: string
    doc: The basename of the index files to write.
    inputBinding:
      position: 2

outputs:
  index_files:
    type: File[]
    outputBinding:
      glob: "*.ebwt"