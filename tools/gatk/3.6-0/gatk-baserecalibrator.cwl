cwlVersion: v1.0
class: CommandLineTool
id: gatk-baserecalibrator
label: "GATK BaseRecalibrator"

baseCommand: [java, -jar, /usr/GenomeAnalysisTK.jar, --analysis_type, BaseRecalibrator]

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
    doc: The output recalibration table file to create
    inputBinding:
      prefix: --out
  knownSites:
    doc: A database of known polymorphic sites
    type:
      type: array
      items: File
      inputBinding:
        prefix: --knownSites
    inputBinding: {}
    secondaryFiles: [ .tbi ]
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
