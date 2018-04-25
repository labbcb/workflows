#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "TrimGalore"
label: "Trim Galore!"
doc: |
  Trim Galore! is a wrapper around Cutadapt and FastQC to consistently apply
  adapter and quality trimming to FastQ files, with extra functionality for RRBS
  data.

baseCommand: [trim_galore]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: welliton/trimgalore:0.4.4

inputs:

# General options

  file_R1:
    type: File
    doc: A single FASTQ file with single-end reads or forward sequences (R1).
    inputBinding:
      position: 1
  file_R2:
    type: File?
    doc: A single FASTQ file with reverse sequences (R2).
    inputBinding:
      position: 2
  quality:
    type: int
    doc: Trim low-quality ends from reads in addition to adapter removal.
    default: 20
    inputBinding:
      prefix: --quality
  phred33:
    type: boolean
    doc: Instructs Cutadapt to use ASCII+33 quality scores as Phred scores.
    default: true
    inputBinding:
      prefix: --phred33
  phred64:
    type: boolean
    doc: Instructs Cutadapt to use ASCII+64 quality scores as Phred scores.
    default: false
    inputBinding:
      prefix: --phred64
  adapter:
    type: string?
    doc: Adapter sequence to be trimmed.
    inputBinding:
      prefix: --adapter
  adapter2:
    type: string?
    doc: |
      Optional adapter sequence to be trimmed off read 2 of paired-end files.
      This option requires --paired to be specified as well.
    inputBinding:
      prefix: --adapter2
  illumina:
    type: boolean
    doc: |
      Adapter sequence to be trimmed is the first 13bp of the Illumina universal
      adapter AGATCGGAAGAGC instead of the default auto-detection of adapter
      sequence.
    default: false
    inputBinding:
      prefix: --illumina
  nextera:
    type: boolean
    doc: |
      Adapter sequence to be trimmed is the first 12bp of the Nextera adapter
      CTGTCTCTTATA instead of the default auto-detection of adapter sequence.
    default: false
    inputBinding:
      prefix: --nextera
  small_rna:
    type: boolean
    doc: |
      Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA
      3' Adapter TGGAATTCTCGG instead of the default auto-detection of adapter
      sequence.
    default: false
    inputBinding:
      prefix: --small_rna
  max_length:
    type: int?
    doc: Discard reads that are longer than bp after trimming.
    inputBinding:
      prefix: --max_length
  stringency:
    type: int
    doc: Overlap with adapter sequence required to trim a sequence.
    default: 1
    inputBinding:
      prefix: --stringency
  error_rate:
    type: float
    doc: |
      Maximum allowed error rate (no. of errors divided by the length of the
      matching region).
    default: 0.1
    inputBinding:
      prefix: -e
  gzip:
    type: boolean
    doc: Compress the output file with gzip.
    default: false
    inputBinding:
      prefix: --gzip
  dont_gzip:
    type: boolean
    doc: Output files won't be compressed with gzip. This overrides --gzip.
    default: false
    inputBinding:
      prefix: --dont_gzip
  length:
    type: int
    doc: |
      Discard reads that became shorter than length INT because of either
      quality or adapter trimming. A value of 0 effectively disables this
      behaviour.
    default: 20
    inputBinding:
      prefix: --length
  max_n:
    type: int?
    doc: |
      The total number of Ns a read may contain before it will be removed
      altogether.
    inputBinding:
      prefix: --max_n
  trim_n:
    type: boolean
    doc: Removes Ns from either side of the read.
    default: false
    inputBinding:
      prefix: --trim-n
  clip_R1:
    type: int?
    doc: Instructs Trim Galore to remove bp from the 5' end of read 1.
    inputBinding:
      prefix: --clip_R1
  clip_R2:
    type: int?
    doc: Instructs Trim Galore to remove bp from the 5' end of read 2.
    inputBinding:
      prefix: --clip_R2
  three_prime_clip_R1:
    type: int?
    doc: |
      Instructs Trim Galore to remove <int> bp from the 3' end of read 1
      (or single-end reads) AFTER adapter/quality trimming has been performed.
    inputBinding:
      prefix: --three_prime_clip_R1
  three_prime_clip_R2:
    type: int?
    doc: |
      Instructs Trim Galore to re move <int> bp from the 3' end of read 2 AFTER
      adapter/quality trimming has been performed.
    inputBinding:
      prefix: --three_prime_clip_R2

# RRBS-specific options (MspI digested material)

  rrbs:
    type: boolean
    doc: Specifies that the input file was an MspI digested RRBS sample.
    default: false
    inputBinding:
      prefix: --rrbs
  non_directional:
    type: boolean
    doc: |
      Selecting this option for non-directional RRBS libraries will screen
      quality-trimmed sequences for CAA or CGA at the start of the read and,
      if found, removes the first two base pairs.
    default: false
    inputBinding:
      prefix: --non_directional
  keep:
    type: boolean
    doc: |
      Keep the quality trimmed intermediate file. If not specified the temporary
      file will be deleted after adapter trimming.
    default: false
    inputBinding:
      prefix: --keep

# Paired-end specific options

  paired:
    type: boolean
    doc: |
      This option performs length trimming of quality/adapter/RRBS trimmed reads
      for paired-end files.
    default: false
    inputBinding:
      prefix: --paired
  trim1:
    type: boolean
    doc: Trims 1 bp off every read from its 3' end.
    default: false
    inputBinding:
      prefix: --trim1

outputs:
  trim_file_R1:
    type: File
    outputBinding:
      glob: $(inputs.file_R1.basename.replace(/\.(fastq|fq)(\.gz)?$/, "_*"))
  trim_file_R2:
    type: File?
    outputBinding:
      glob: ["*_val_2.fq", "*_val_2.fq.gz"]
  stats_file_R1:
    type: File
    outputBinding:
      glob: "$(inputs.file_R1.basename)_trimming_report.txt"
  stats_file_R2:
    type: File?
    outputBinding:
      glob: |
        ${
          return inputs.file_R2 ? inputs.file_R2.basename.concat("_trimming_report.txt") : "";
        }
