cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:

# General
  files: File[]
  workers: int
  groups: string[]

# Quality assessment of raw data
  reportFileRaw: string
  rdsFileRaw: string

# Trim raw sequences
  small_rna: boolean
  length: int
  max_length: int
  max_n: int
  dont_gzip: boolean

# Quality assessment of trimmed data
  reportFileTrim: string
  rdsFileTrim: string

# Build genome index
  genome_files: File[]
  ebwt_base: string

# Align trimmed sequences
  seedmms: int
  all: boolean
  best: boolean
  strata: boolean
  sam: boolean

# Count aligned sequences
  count_file_name: string
  annotation: File
  feature_type: string
  attribute_type: string
  feature_level: boolean
  overlaps: boolean
  multi_mapping: boolean
  mapping_quality: int

outputs:
  count_file:
    type: File
    outputSource: count/output

# Reports and summaries
  qa_report_raw:
    type: File
    outputSource: qa_raw/qc_report
  trim_stats_files:
    type: File[]
    outputSource: trim/stats_file_R1
  qa_report_trim:
    type: File
    outputSource: qa_trim/qc_report
  count_summary:
    type: File
    outputSource: count/summary

# Intermediate results
#  qa_rds_raw:
#    type: File
#    outputSource: qa_raw/qc_rds
#  trim_files:
#    type: File[]
#    outputSource: trim/trim_file_R1
#  qa_rds_trim:
#    type: File
#    outputSource: qa_trim/qc_rds
#  align_files:
#    type: File[]
#    outputSource: align/output

steps:
  qa_raw:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Frqc/versions/v1.12.0/plain-CWL/descriptor
    in:
      files: files
      groups: groups
      workers: workers
      reportFile: reportFileRaw
      rdsFile: rdsFileRaw
    out: [qc_report, qc_rds]
  trim:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Ftrimgalore/versions/v0.4.4/plain-CWL/descriptor
    scatter: file_R1
    in:
      file_R1: files
      small_rna: small_rna
      length: length
      max_length: max_length
      max_n: max_n
      dont_gzip: dont_gzip
    out: [trim_file_R1, stats_file_R1]
  qa_trim:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Frqc/versions/v1.12.0/plain-CWL/descriptor
    in:
      files: trim/trim_file_R1
      groups: groups
      workers: workers
      reportFile: reportFileTrim
      rdsFile: rdsFileTrim
    out: [qc_report, qc_rds]
  build_index:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Fbowtie%2Fbuild/versions/v1.2.0/plain-CWL/descriptor
    in:
      genome_files: genome_files
      ebwt_base: ebwt_base
    out: [index_files]
  align:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Fbowtie/versions/v1.2.0/plain-CWL/descriptor
    scatter: [file]
    in:
      index_files: build_index/index_files
      ebwt: ebwt_base
      file: trim/trim_file_R1
      seedmms: seedmms
      all: all
      best: best
      strata: strata
      sam: sam
      threads: workers
      hit:
        valueFrom: $(inputs.file.basename).sam
    out: [output]
  count:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Fsubread%2Ffeaturecounts/versions/v1.5.2/plain-CWL/descriptor
    in:
      files: align/output
      annotation: annotation
      feature_type: feature_type
      attribute_type: attribute_type
      feature_level: feature_level
      overlaps: overlaps
      multi_mapping: multi_mapping
      mapping_quality: mapping_quality
      nthreads: workers
      output_file: count_file_name
    out: [output, summary]
