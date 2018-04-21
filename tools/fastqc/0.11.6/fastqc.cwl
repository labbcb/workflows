#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "fastqc" 
label: "A high throughput sequence QC analysis tool"
doc: |
  FastQC aims to provide a simple way to do some quality control checks on raw
  sequence data coming from high throughput sequencing pipelines. It provides a
  modular set of analyses which you can use to give a quick impression of
  whether your data has any problems of which you should be aware before doing
  any further analysis. 

baseCommand: fastqc

requirements:
  InitialWorkDirRequirement:
    listing: [$(inputs.file)]
  DockerRequirement:
    dockerPull: welliton/fastqc:v0.11.6

inputs:

  file:
    type: File
    doc: Input file.
    inputBinding:
      position: 2
  casava:
    type: boolean
    default: false
    inputBinding:
      prefix: --casava
      position: 1
    doc: Files come from raw casava output.
  nano:
    type: boolean
    default: false
    inputBinding:
      prefix: --nano
      position: 1
    doc: Files come from naopore sequences and are in fast5 format.
  nofilter:
    type: boolean
    default: false
    inputBinding:
      prefix: --nofilter
      position: 1
    doc: |
      If running with --casava then don't remove read flagged by casava as poor
      quality when performing the QC analysis.
  extract:
    type: boolean
    default: false
    inputBinding:
      prefix: --extract
      position: 1
    doc: |
      If set then the zipped output file will be uncompressed in the same 
      directory after it has been created.
  nogroup:
    type: boolean
    default: false
    inputBinding:
      prefix: --nogroup
      position: 1
    doc: Disable grouping of bases for reads >50bp.
  min_length:
    type: int?
    inputBinding:
      prefix: --min_length
      position: 1
    doc: |
      Sets an artificial lower limit on the length of the sequence to be shown
      in the report.
  format:
    type: string?
    inputBinding:
      prefix: --format
      position: 1
    doc: |
      Bypasses the normal sequence file format detection and forces the program
      to use the specified format. Valid formats are bam, sam, bam_mapped,
      sam_mapped and fastq.
  threads:
    type: int
    default: 1
    inputBinding:
      prefix: --threads
      position: 1
    doc: Specifies the number of files which can be processed simultaneously.
  contaminants:
    type: File?
    inputBinding:
      prefix: --contaminants
      position: 1
    doc: |
      Specifies a non-default file which contains the list of contaminants to
      screen overrepresented sequences against.
  adapters:
    type: File?
    inputBinding:
      prefix: --adapters
      position: 1
    doc: |
      Specifies a non-default file which contains the list of adapter sequences
      which will be explicity searched against the library.
  limits:
    type: File?
    inputBinding:
      prefix: --limits
      position: 1
    doc: |
      Specifies a non-default file which contains a set of criteria which will
      be used to determine the warn/error limits for the various modules.
  kmers:
    type: int
    default: 7
    inputBinding:
      prefix: --kmers
      position: 1
    doc: |
      Specifies the length of Kmer to look for in the Kmer content module.
      Specified Kmer length must be between 2 and 10. Default length is 7 if not
      specified.
  quiet:
    type: boolean
    default: true
    inputBinding:
      prefix: --quiet
      position: 1
    doc: Supress all progress messages on stdout and only report errors.

outputs:
  report:
    type: File
    outputBinding:
      glob: "*_fastqc.html"
  report_zip:
    type: File
    outputBinding:
      glob: "*_fastqc.zip"