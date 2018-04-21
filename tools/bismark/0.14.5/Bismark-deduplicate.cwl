cwlVersion: v1.0
class: CommandLineTool

id: "bismark_deduplicate"
label: "Bismark deduplicate"
doc: |
  Bismark is a program to map bisulfite treated sequencing reads to a genome of
  interest and perform methylation calls in a single step.

dct:creator:
  "@id": "https://orcid.org/0000-0001-8377-1836"
  foaf:name: Welliton Souza
  foaf:mbox: "mailto:well309@gmail.com"

baseCommand: [deduplicate_bismark]

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.file)
  DockerRequirement:
    dockerPull: welliton/bismark:v0.14.5

inputs:
  file:
    type: File
    inputBinding:
      valueFrom: $(self.basename)
  single:
    type: boolean?
    doc: Deduplicate single-end Bismark files (default format SAM).
    inputBinding:
      prefix: --single
  paired:
    type: boolean?
    doc: Deduplicate paired-end Bismark files (default format SAM).
    inputBinding:
      prefix: --paired
  vanilla:
    type: boolean
    doc: The input file is in the old custom Bismark format and not in SAM format.
    default: false
    inputBinding:
      prefix: --vanilla
  barcode:
    type: string?
    doc: |
      In addition to chromosome, start position and orientation this will also
      take a potential barcode into consideration while deduplicating. The
      barcode needs to be the last element of the read ID and separated by a
      ':', e.g.: MISEQ:14:000000000-A55D0:1:1101:18024:2858_1:N:0:CTCCT.
    inputBinding:
      prefix: --barcode
  bam:
    type: boolean
    doc: |
      The output will be written out in BAM format instead of the default SAM
      format. This script will attempt to use the path to Samtools that was
      specified with '--samtools_path', or, if it hasn't been specified, attempt
      to find Samtools in the PATH. If no installation of Samtools can be found,
      the SAM output will be compressed with GZIP instead (yielding a .sam.gz
      output file).
    default: false
    inputBinding:
      prefix: --bam
  samtools_path:
    type: string?
    doc: |
      The path to your Samtools installation, e.g. /home/user/samtools/.
      Does not need to be specified explicitly if Samtools is in the PATH
      already.
    inputBinding:
      prefix: --samtools_path

outputs:
  output:
    type: File
    outputBinding:
      glob: "*deduplicated*"
  report:
    type: File
    outputBinding:
      glob: "*_report.txt"
