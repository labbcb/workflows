cwlVersion: v1.0
class: CommandLineTool
id: gatk-printreads
label: "GATK PrintReads"

baseCommand: [java, -jar, /usr/GenomeAnalysisTK.jar, --analysis_type, PrintReads]

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference_sequence)
      - $(inputs.reference_index)
      - $(inputs.reference_dictionary)
  DockerRequirement:
    dockerPull: broadinstitute/gatk3:3.6-0

inputs:
  input_file:
    type: File
    doc: Input file containing sequence data (BAM or CRAM)
    inputBinding:
      prefix: --input_file
    secondaryFiles: [ ^.bai ]
  reference_sequence:
    type: File
    doc: Reference sequence file
    inputBinding:
      prefix: --reference_sequence
  reference_index:
    type: File
  reference_dictionary:
    type: File
  out:
    type: string
    doc: Write output to this BAM filename
    inputBinding:
      prefix: --out
  intervals:
    type: File?
    doc: One or more genomic intervals over which to operate
    inputBinding:
      prefix: --intervals
  interval_padding:
    type: int?
    doc: Amount of padding (in bp) to add to each interval
    inputBinding:
      prefix: --interval_padding
  BQSR:
    type: File?
    doc: |
      Input covariates table file for on-the-fly base quality score
      recalibration.
    inputBinding:
      prefix: --BQSR

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)
    secondaryFiles:
      - ^.bai
