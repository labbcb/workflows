cwlVersion: v1.0
class: CommandLineTool
id: picard-buildbamindex
label: "Picard BuildBamIndex"
doc: |
  Generates a BAM index ".bai" file. This tool creates an index file for the
  input BAM that allows fast look-up of data in a BAM file, lke an index on a
  database. Note that this tool cannot be run on SAM files, and that the input
  BAM file must be sorted in coordinate order.

baseCommand: [java, -jar, /usr/picard.jar, BuildBamIndex]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.input) ]

hints:
  DockerRequirement:
    dockerPull: welliton/picard:2.6.0

inputs:
  input:
    type: File
    doc: |
      A BAM file or GA4GH URL to process. Must be sorted in coordinate order.
    inputBinding:
      prefix: INPUT=
      separate: false

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.input.basename)
    secondaryFiles:
      - ^.bai
