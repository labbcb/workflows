cwlVersion: v1.0
class: Workflow

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:

# General
  files_R1: File[]
  files_R2: File[]
  workers: int
  groups_R1: string[]
  groups_R2: string[]
  paired: boolean

# Quality assessment of raw data
  sample: boolean
  reportFileRaw: string
  rdsFileRaw: string

# Trim raw sequences
  trim1: boolean
  illumina: boolean

# Quality assessment of trimmed data
  reportFileTrim: string
  rdsFileTrim: string

# Prepare reference genome
  genome_files: File[]

# Align trimmed sequences
  seedmms: int
  bowtie1: boolean
  unmapped: boolean
  ambigous: boolean

# Deduplicate aligned sequences
  bam: boolean

# Methylation extraction
  bedGraph: boolean

outputs:
  coverage_files:
    type: File[]?
    outputSource: methylation-extraction/coverage
  bedGraph_files:
    type: File[]?
    outputSource: methylation-extraction/bedGraph_file

# Reports and summaries
  qa_report_raw:
    type: File
    outputSource: qa_raw/qc_report
  trim_stats_R1:
    type: File[]
    outputSource: trim/stats_file_R1
  trim_stats_R2:
    type: File[]
    outputSource: trim/stats_file_R2
  qa_report_trim:
    type: File
    outputSource: qa_trim/qc_report
  align_report:
    type: File[]
    outputSource: align/report
  deduplicate_report:
    type: File[]
    outputSource: deduplicate/report
  mbias_plot_R1:
    type: File[]?
    outputSource: methylation-extraction/mbias_plot_R1
  mbias_plot_R2:
    type: File[]?
    outputSource: methylation-extraction/mbias_plot_R2
  mbias_report:
    type: File[]?
    outputSource: methylation-extraction/mbias_report
  methylation_extration_report:
    type: File[]
    outputSource: methylation-extraction/report_file
#  qa_rds_raw:
#    type: File
#    outputSource: qa_raw/qc_rds
#  qa_rds_trim:
#    type: File
#    outputSource: qa_trim/qc_rds

# Intermediate results
#  trim_files:
#    type: File[]
#    outputSource: [trim/trim_file_R1, trim/trim_file_R2]
#    linkMerge: merge_flattened
#  align_files:
#    type: File[]
#    outputSource: align/output
#  deduplicate_files:
#    type: File[]
#    outputSource: deduplicate/output


steps:
  qa_raw:
    run: https://raw.githubusercontent.com/labbcb/docker-rqc/master/Rqc.cwl
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
    out: [qc_report, qc_rds]
  trim:
    run: https://raw.githubusercontent.com/labbcb/tool-trimgalore/master/TrimGalore.cwl
    scatter: [file_R1, file_R2]
    scatterMethod: dotproduct
    in:
      file_R1: files_R1
      file_R2: files_R2
      paired: paired
      trim1: trim1
      illumina: illumina
    out: [trim_file_R1, trim_file_R2, stats_file_R1, stats_file_R2]
  qa_trim:
    run: https://raw.githubusercontent.com/labbcb/docker-rqc/master/Rqc.cwl
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
    out: [qc_report, qc_rds]
  prepare-genome:
    run: https://raw.githubusercontent.com/labbcb/tool-bismark/master/Bismark-genome-preparation.cwl
    in:
      genome_files: genome_files
      bowtie1: bowtie1
    out: [genome_folder]
  align:
    run: https://raw.githubusercontent.com/labbcb/tool-bismark/master/Bismark.cwl
    scatter: [mate1, mate2]
    scatterMethod: dotproduct
    in:
      genome_files: genome_files
      genome_folder: prepare-genome/genome_folder
      mate1: trim/trim_file_R1
      mate2: trim/trim_file_R2
      seedmms: seedmms
      bowtie1: bowtie1
      unmapped: unmapped
      ambigous: ambigous
      multicore: workers
    out: [output, report]
  deduplicate:
    run: https://raw.githubusercontent.com/labbcb/tool-bismark/master/Bismark-deduplicate.cwl
    scatter: [file]
    in:
      file: align/output
      paired: paired
      bam: bam
    out: [output, report]
  methylation-extraction:
    run: https://raw.githubusercontent.com/labbcb/tool-bismark/master/Bismark-methylation-extractor.cwl
    scatter: [file]
    in:
      file: deduplicate/output
      paired_end: paired
      bedGraph: bedGraph
      multicore: workers
    out: [bedGraph_file, coverage, mbias_plot_R1, mbias_plot_R2, mbias_report, report_file]