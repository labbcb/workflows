cwlVersion: v1.0
class: CommandLineTool
id: gatk-haplotypecaller
label: "GATK HaplotypeCaller"

baseCommand: [java, -jar, /usr/GenomeAnalysisTK.jar, --analysis_type, HaplotypeCaller]

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
    doc: File to which variants should be written
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
  genotyping_mode:
    type: string?
    doc: Specifies how to determine the alternate alleles to use for genotyping
    inputBinding:
      prefix: --genotyping_mode
  stand_call_conf:
    type: int?
    doc: |
      The minimum phred-scaled confidence threshold at which variants should be
      called
    inputBinding:
      prefix: -stand_call_conf
  stand_emit_conf:
    type: int?
    doc: |
      The minimum phred-scaled confidence threshold at which variants should be
      emitted (and filtered with LowQual if less than the calling threshold)
    inputBinding:
      prefix: -stand_emit_conf
  emitRefConfidence:
    type: string?
    doc: Mode for emitting reference confidence scores
    inputBinding:
      prefix: --emitRefConfidence
  variant_index_type:
    type: string?
    doc: Type of IndexCreator to use for VCF/BCF indices
    inputBinding:
      prefix: --variant_index_type
  variant_index_parameter:
    type: int?
    doc: Parameter to pass to the VCF/BCF IndexCreator
    inputBinding:
      prefix: --variant_index_parameter

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)
    secondaryFiles: [ .idx ]
