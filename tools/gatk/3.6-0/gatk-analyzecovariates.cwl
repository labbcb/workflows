cwlVersion: v1.0
class: CommandLineTool
id: gatk-analyzecovariates
label: "GATK AnalyzeCovariates"

baseCommand: [java, -jar, /usr/GenomeAnalysisTK.jar, --analysis_type, AnalyzeCovariates]

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference_sequence)
      - $(inputs.reference_index)
      - $(inputs.reference_dictionary)
  DockerRequirement:
    dockerPull: broadinstitute/gatk3:3.6-0

inputs:
  beforeReportFile:
    type: File?
    doc: File containing the recalibration tables from the first pass.
    inputBinding:
      prefix: --beforeReportFile
  afterReportFile:
    type: File?
    doc: File containing the recalibration tables from the second pass.
    inputBinding:
      prefix: --afterReportFile
  reference_sequence:
    type: File
    doc: Reference sequence file
    inputBinding:
      prefix: --reference_sequence
  reference_index:
    type: File
  reference_dictionary:
    type: File
  plotsReportFile:
    type: string
    doc: Output report file name.
    inputBinding:
      prefix: --plotsReportFile
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

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.plotsReportFile)
