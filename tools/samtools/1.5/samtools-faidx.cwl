cwlVersion: v1.0
class: CommandLineTool
id: samtools-faidx
label: "Samtools faidx"
doc: |
  Index reference sequence in the FASTA format or extract subsequence from
  indexed reference sequence. If no region is specified, faidx will index the
  file and create <ref.fasta>.fai on the disk. If regions are specified, the
  subsequences will be retrieved and printed to stdout in the FASTA format.
  The input file can be compressed in the BGZF format.
  The sequences in the input file should all have different names. If they do not,
  indexing will emit a warning about duplicate sequences and retrieval will only
  produce subsequences from the first sequence with the duplicated name.

baseCommand: [samtools, faidx]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.ref_fasta) ]
  DockerRequirement:
    dockerPull: welliton/samtools:1.5

inputs:
  ref_fasta:
    type: File
    doc: Reference genome file in FASTA format.
    inputBinding: {}

outputs:
  index_file:
    type: File
    outputBinding:
      glob: $(inputs.ref_fasta.basename).fai
