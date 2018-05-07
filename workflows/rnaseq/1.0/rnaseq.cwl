cwlVersion: v1.0
class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}

inputs:

# General
  files_R1: File[]
  files_R2: File[]
  workers: int
  gtf_file: File

# Quality assessment of raw data
  groups_R1: string[]
  groups_R2: string[]
  sample: boolean
  reportFileRaw: string
  rdsFileRaw: string

# Trim raw sequencing reads
  illumina: boolean
  paired: boolean
  dont_gzip: boolean

# QA trim
  reportFileTrim: string
  rdsFileTrim: string

# Align
  genomeFastaFiles: File[]
  outSAMtype: string[]?
  outFileNamePrefix: string[]

# Count
  f: string
  stranded: string

outputs:
# Quality control
  qa_report_raw:
    type: File
    outputSource: rqc_raw/qc_report
  fastqc_raw_reports:
    type: File[]
    outputSource: fastqc_raw/report
  trim_stats_files_R1:
    type: File[]
    outputSource: trim/stats_file_R1
  trim_stats_files_R2:
    type: File[]
    outputSource: trim/stats_file_R2
  fastqc_raw_reports_zip:
    type: File[]
    outputSource: fastqc_raw/report_zip
  qa_report_trim:
    type: File
    outputSource: rqc_trim/qc_report
  fastqc_trim_reports:
    type: File[]
    outputSource: fastqc_trim/report
  fastqc_trim_reports_zip:
    type: File[]
    outputSource: fastqc_trim/report_zip

# Count
  count_files:
    type: File[]
    outputSource: count/output

# Temporary files
  trim_files_R1:
    type: File[]
    outputSource: trim/trim_file_R1
  trim_files_R2:
    type: File[]
    outputSource: trim/trim_file_R2
  index_files:
    type: File[]
    outputSource: build_index/index_files
  align_files:
    type: File[]
    outputSource: align/output

steps:
  rqc_raw:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Frqc/versions/v1.12.0/plain-CWL/descriptor
    in:
      files:
        source: [files_R1, files_R2]
        linkMerge: merge_flattened
      groups:
        source: [groups_R1, groups_R2]
        linkMerge: merge_flattened
      sample: sample
      workers: workers
      reportFile: reportFileRaw
      rdsFile: rdsFileRaw
    out: [qc_report]
  fastqc_raw:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Ffastqc/versions/v0.11.6/plain-CWL/descriptor
    scatter: [file]
    in:
      file:
        source: [files_R1, files_R2]
        linkMerge: merge_flattened
    out: [report, report_zip]
  trim:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Ftrimgalore/versions/v0.4.4/plain-CWL/descriptor
    scatter: [file_R1, file_R2]
    scatterMethod: dotproduct
    in:
      file_R1: files_R1
      file_R2: files_R2
      illumina: illumina
      paired: paired
      dont_gzip: dont_gzip
    out: [trim_file_R1, stats_file_R1, trim_file_R2, stats_file_R2]
  rqc_trim:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Frqc/versions/v1.12.0/plain-CWL/descriptor
    in:
      files:
        source: [trim/trim_file_R1, trim/trim_file_R2]
        linkMerge: merge_flattened
      groups:
        source: [groups_R1, groups_R2]
        linkMerge: merge_flattened
      sample: sample
      workers: workers
      reportFile: reportFileTrim
      rdsFile: rdsFileTrim
    out: [qc_report]
  fastqc_trim:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Ffastqc/versions/v0.11.6/plain-CWL/descriptor
    scatter: [file]
    in:
      file:
        source: [trim/trim_file_R1, trim/trim_file_R2]
        linkMerge: merge_flattened
    out: [report, report_zip]
  build_index:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Fstar%2Fgenomegenerate/versions/v2.5.3a/plain-CWL/descriptor
    in:
      genomeFastaFiles: genomeFastaFiles
      sjdbGTFfile: gtf_file
    out: [index_files]
  align:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Fstar%2Falignreads/versions/v2.5.3a/plain-CWL/descriptor
    scatter: [readFilesIn_R1, readFilesIn_R2, outFileNamePrefix]
    scatterMethod: dotproduct
    in:
      readFilesIn_R1: trim/trim_file_R1
      readFilesIn_R2: trim/trim_file_R2
      index_files: build_index/index_files
      outSAMtype: outSAMtype
      outFileNamePrefix: outFileNamePrefix
    out: [output]
  count:
    run: https://dockstore.org:8443/api/ga4gh/v1/tools/registry.hub.docker.com%2Fwelliton%2Fhtseq/versions/v0.9.1/plain-CWL/descriptor
    scatter: [input]
    in:
      input: align/output
      gtf_file: gtf_file
      f: f
      stranded: stranded
    out: [output]
