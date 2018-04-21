#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
id: "featureCounts" 
label: "Ultrafast and accurate read summarization program."
doc: |
  featureCounts is a highly efficient general-purpose read summarization program
  that counts mapped reads for genomic features such as genes, exons, promoter,
  gene bodies, genomic bins and chromosomal locations. It can be used to count
  both RNA-seq and genomic DNA-seq reads.
  
baseCommand: featureCounts

requirements:
  DockerRequirement:
    dockerPull: welliton/subread:v1.5.2

inputs:

# Mandatory arguments

  annotation:
    type: File
    inputBinding:
      prefix: -a
    doc: Name of an annotation file.
  output_file:
    type: string
    inputBinding:
      prefix: -o
    doc: Name of the output file including read counts.
  files:
    type: File[]
    inputBinding: {}
    doc: A list of SAM or BAM format files.

# Annotation

  format:
    type: string
    default: GTF
    inputBinding:
      prefix: -F
    doc: Specify format of the provided annotation file.
  feature_type:
    type: string
    default: exon
    inputBinding:
      prefix: -t
    doc: Specify feature type in GTF annotation.
  attribute_type:
    type: string
    default: gene_id
    inputBinding:
      prefix: -g
    doc: Specify attribute type in GTF annotation.

# Level of summarization

  feature_level:
    type: boolean
    default: false
    inputBinding:
      prefix: -f
    doc: Perform read counting at feature level.
  
# Overlap between reads and features

  overlaps:
    type: boolean
    default: false
    inputBinding:
      prefix: -O
    doc: |
      Assign reads to all their overlapping meta-features
      (or features if feature_level is true).

# Multi-mapping reads

  multi_mapping:
    type: boolean
    default: false
    inputBinding:
      prefix: -M
    doc: Multi-mapping reads will also be counted.

# Read filtering

  mapping_quality:
    type: int
    default: 0
    inputBinding:
      prefix: -Q
    doc: |
      The minimum mapping quality score a read must satisfy in order to be
      counted.

# Number of CPU threads

  nthreads:
    type: int
    default: 1
    inputBinding:
      prefix: -T
    doc: Number of the threads. 1 by default.

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_file)
  summary:
    type: File
    outputBinding:
      glob: "*.summary"
