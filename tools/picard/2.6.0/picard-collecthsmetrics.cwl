cwlVersion: v1.0
class: CommandLineTool
id: picard-collecthsmetrics
label: "Picard CollectHsMetrics"
doc: |
  Collects hybrid-selection (HS) metrics for a SAM or BAM file. This tool takes
  a SAM/BAM file input and collects metrics that are specific for sequence
  datasets generated through hybrid-selection. Hybrid-selection (HS) is the most
  commonly used technique to capture exon-specific sequences for targeted
  sequencing experiments such as exome sequencing; for more information, please
  see the corresponding GATK Dictionary entry
  (http://www.broadinstitute.org/gatk/guide/article?id=6331).
  This tool requires an aligned SAM or BAM file as well as bait and target
  interval files in Picard interval_list format. You should use the bait and
  interval files that correspond to the capture kit that was used to generate
  the capture libraries for sequencing, which can generally be obtained from the
  kit manufacturer. If the baits and target intervals are provided in BED format,
  you can convert them to the Picard interval_list format using Picards
  BedToInterval tool.
  If a reference sequence is provided, this program will calculate both
  AT_DROPOUT and GC_DROPOUT metrics. Dropout metrics are an attempt to measure
  the reduced representation of reads, in regions that deviate from 50% G/C
  content. This reduction in the number of aligned reads is due to the increased
  numbers of errors associated with sequencing regions with excessive or
  deficient numbers of G/C bases, ultimately leading to poor mapping
  efficiencies and lowcoverage in the affected regions.
  If you are interested in getting G/C content and mean sequence depth
  information for every target interval, use the PER_TARGET_COVERAGE option.
  Metrics labeled as percentages are actually expressed as fractions!

baseCommand: [java, -jar, /usr/picard.jar, CollectHsMetrics]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.reference_sequence), $(inputs.reference_index) ]

hints:
  DockerRequirement:
    dockerPull: welliton/picard:2.6.0

inputs:
  input:
    type: File
    doc: An aligned SAM or BAM file.
    inputBinding:
      prefix: INPUT=
      separate: false
  output:
    type: string
    doc: The output file to write the metrics to.
    inputBinding:
      prefix: OUTPUT=
      separate: false
  target_intervals:
    type: File
    doc: |
      An interval list file that contains the locations of the targets. Default
      value null. This option must be specified at least 1 times.
    inputBinding:
      prefix: TARGET_INTERVALS=
      separate: false
  bait_intervals:
    type: File
    doc: |
      An interval list file that contains the locations of the baits used.
      Default value null. This option must be specified at least 1 times.
    inputBinding:
      prefix: BAIT_INTERVALS=
      separate: false
  reference_sequence:
    type: File?
    inputBinding:
      prefix: REFERENCE_SEQUENCE=
      separate: false
  reference_index:
    type: File?
  bait_set_name:
    type: string?
    doc: |
      Bait set name. If not provided it is inferred from the filename of the
      bait intervals. Default value null.
    inputBinding:
      prefix: BAIT_SET_NAME=
      separate: false
  metric_accumulation_level:
    type: string?
    doc: |
      The level(s) at which to accumulate metrics. Default value [ALL_READS].
      This option can be set to 'null' to clear the default value. Possible
      values {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be
      specified 0 or more times. This option can be set to 'null' to clear the
      default list.
    inputBinding:
      prefix: METRIC_ACCUMULATION_LEVEL=
      separate: false
  per_target_coverage:
    type: string?
    doc: |
      An optional file to output per target coverage information to.
      Default value null.
    inputBinding:
      prefix: PER_TARGET_COVERAGE=
      separate: false
  per_base_coverage:
    type: string?
    doc: |
      An optional file to output per base coverage information to. The per-base
      file contains one line per target base and can grow very large. It is not
      recommended for use with large target sets. Default value null.
    inputBinding:
      prefix: PER_BASE_COVERAGE=
      separate: false
  near_distance:
    type: int?
    doc: |
      The maximum distance between a read and the nearest probe/bait/amplicon
      for the read to be considered 'near probe' and included in percent
      selected. Default value 250. This option can be set to 'null' to clear
      the default value.
    inputBinding:
      prefix: NEAR_DISTANCE=
      separate: false
  minimum_mapping_quality:
    type: int?
    doc: |
      Minimum mapping quality for a read to contribute coverage. Default value
      20. This option can be set to 'null' to clear the default value.
    inputBinding:
      prefix: MINIMUM_MAPPING_QUALITY=
      separate: false
  minimum_basequality:
    type: int?
    doc: |
      Minimum base quality for a base to contribute coverage. Default value 20.
      This option can be set to 'null' to clear the default value.
    inputBinding:
      prefix: MINIMUM_BASE_QUALITY=
      separate: false
  clip_overlapping_reads:
    type: string?
    doc: |
      True if we are to clip overlapping reads, false otherwise. Default value
      true. This option can be set to 'null' to clear the default value.
      Possible values {true, false}
    inputBinding:
      prefix: CLIP_OVERLAPPING_READS=
      separate: false
  coverage_cap:
    type: int?
    doc: |
       Parameter to set a max coverage limit for Theoretical Sensitivity
       calculations. Default is 200. This option can be set to 'null' to clear
       the default value.
    inputBinding:
      prefix: COVERAGE_CAP=
      separate: false
  sample_size:
    type: int?
    doc: |
      Sample Size used for Theoretical Het Sensitivity sampling. Default is
      10000. This option can be set to 'null' to clear the default value.
    inputBinding:
      prefix: SAMPLE_SIZE=
      separate: false
  verbosity:
    type: string
    default: INFO
    doc: |
      Control verbosity of logging. Default value INFO. This option can be set
      to 'null' to clear the default value.
      Possible values {ERROR, WARNING, INFO, DEBUG}
    inputBinding:
      prefix: VERBOSITY=
      separate: false

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output)
  per_target_coverage_file:
    type: File?
    outputBinding:
      glob: $(inputs.per_target_coverage)
