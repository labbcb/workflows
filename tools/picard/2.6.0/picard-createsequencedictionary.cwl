cwlVersion: v1.0
class: CommandLineTool
id: picard-createsequencedictionary
label: "Picard CreateSequenceDictionary"

baseCommand: CreateSequenceDictionary

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.reference) ]

hints:
  DockerRequirement:
    dockerPull: welliton/picard:2.6.0

inputs:
  reference:
    type: File
    doc: Input reference fasta or fasta.gz
    inputBinding:
      prefix: REFERENCE=
      separate: false
  output:
    type: string
    doc: |
       Output SAM or BAM file containing only the sequence dictionary.
    inputBinding:
      prefix: OUTPUT=
      separate: false
  genome_assembly:
    type: string?
    doc: Put into AS field of sequence dictionary entry if supplied.
    inputBinding:
      prefix: GENOME_ASSEMBLY=
      separate: false
  uri:
    type: string?
    doc: |
      Put into UR field of sequence dictionary entry. If not supplied, input
      reference file is used.
    inputBinding:
      prefix: URI=
      separate: false
  species:
    type: string?
    doc: Put into SP field of sequence dictionary entry
    inputBinding:
      prefix: SPECIES=
      separate: false
  truncate_names_at_whitespace:
    type: string?
    doc: |
      Make sequence name the first word from the > line in the fasta file. By
      default the entire contents of the > line is used, excluding leading and
      trailing whitespace.  Default value true. This option can be set to 'null'
      to clear the default value. Possible values {true, false}.
    inputBinding:
      prefix: TRUNCATE_NAMES_AT_WHITESPACE=
      separate: false
  num_sequences:
    type: int?
    doc: |
      Stop after writing this many sequences. For testing. Default value
      2147483647. This option can be set to 'null' to clear the default value.
    inputBinding:
      prefix: NUM_SEQUENCES=
      separate: false

outputs:
  output_file:
    type: File
    outputBinding:
      glob: "*.dict"
