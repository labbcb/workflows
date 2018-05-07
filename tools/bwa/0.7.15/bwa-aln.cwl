cwlVersion: v1.0
class: CommandLineTool
id: bwa-aln
label: "BWA aln"
doc: Gapped/ungapped alignment.

baseCommand: [bwa, aln]

hints:
  DockerRequirement:
    dockerPull: welliton/bwa:0.7.15

inputs:
  prefix:
    type: Directory
    inputBinding:
      position: 1
  file:
    type: File
    inputBinding:
      position: 2
  n:
    type: float
    default: 0.04
    doc: Max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
    inputBinding:
      prefix: -n
  o:
    type: int
    default: 1
    doc: Maximum number or fraction of gap opens [1]
    inputBinding:
      prefix: -o
  e:
    type: int
    default: -1
    doc: Maximum number of gap extensions, -1 for disabling long gaps [-1]
    inputBinding:
      prefix: -e
  i:
    type: int
    default: 5
    doc: Do not put an indel within INT bp towards the ends [5]
    inputBinding:
      prefix: -i
  d:
    type: int
    default: 10
    doc: Maximum occurrences for extending a long deletion [10]
    inputBinding:
      prefix: -d
  l:
    type: int
    default: 32
    doc: Seed length [32]
    inputBinding:
      prefix: -l
  k:
    type: int
    default: 2
    doc: Maximum differences in the seed [2]
    inputBinding:
      prefix: -k
  m:
    type: int
    default: 2000000
    doc: Maximum entries in the queue [2000000]
    inputBinding:
      prefix: -m
  t:
    type: int
    default: 1
    doc: Number of threads [1]
    inputBinding:
      prefix: -t
  M:
    type: int
    default: 3
    doc: Mismatch penalty [3]
    inputBinding:
      prefix: -M
  O:
    type: int
    default: 11
    doc: Gap open penalty [11]
    inputBinding:
      prefix: -O
  E:
    type: int
    default: 4
    doc: Gap extension penalty [4]
    inputBinding:
      prefix: -E
  R:
    type: int
    default: 30
    doc: Stop searching when there are >INT equally best hits [30]
    inputBinding:
      prefix: -R
  q:
    type: int
    default: 0
    doc: Quality threshold for read trimming down to 35bp [0]
    inputBinding:
      prefix: -q
  f:
    type: string?
    doc: File to write output to instead of stdout
    inputBinding:
      prefix: -f
  B:
    type: int?
    doc: Length of barcode
    inputBinding:
      prefix: -B
  L:
    type: boolean
    default: false
    doc: Log-scaled gap penalty for long deletions
    inputBinding:
      prefix: -L
  N:
    type: boolean
    default: false
    doc: Non-iterative mode search for all n-difference hits (slooow)
    inputBinding:
      prefix: -N
  I:
    type: boolean
    default: false
    doc: The input is in the Illumina 1.3+ FASTQ-like format
    inputBinding:
      prefix: -I
  b:
    type: boolean
    default: false
    doc: The input read file is in the BAM format
    inputBinding:
      prefix: -b
  single-end:
    type: boolean
    default: false
    doc: Use single-end reads only (effective with -b)
    inputBinding:
      prefix: "-0"
  first-read:
    type: boolean
    default: false
    doc: Use the 1st read in a pair (effective with -b)
    inputBinding:
      prefix: "-1"
  second-read:
    type: boolean
    default: false
    doc: Use the 2nd read in a pair (effective with -b)
    inputBinding:
      prefix: "-2"
  casava:
    type: boolean
    default: false
    doc: Filter Casava-filtered sequences
    inputBinding:
      prefix: -Y

outputs:
  output:
    type: stdout
