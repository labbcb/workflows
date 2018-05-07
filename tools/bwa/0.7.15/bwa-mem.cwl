cwlVersion: v1.0
class: CommandLineTool
id: bwa-mem
label: "BWA mem"
doc: BWA-MEM algorithm

baseCommand: [bwa, mem]

hints:
  DockerRequirement:
    dockerPull: welliton/bwa:0.7.15

inputs:
  idxbase:
    type: File
    inputBinding:
      position: 1
  file_R1:
    type: File
    inputBinding:
      position: 2
  file_R2:
    type: File?
    inputBinding:
      position: 3
  output:
    type: string
    default: output.sam

# Algorithm options
  t:
    type: int?
    doc: number of threads [1]
    inputBinding:
      prefix: -t
  k:
    type: int?
    doc: minimum seed length [19]
    inputBinding:
      prefix: -k
  w:
    type: int?
    doc: band width for banded alignment [100]
    inputBinding:
      prefix: -w
  d:
    type: int?
    doc: off-diagonal X-dropoff [100]
    inputBinding:
      prefix: -d
  r:
    type: float?
    doc: look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
    inputBinding:
      prefix: -r
  y:
    type: int?
    doc: seed occurrence for the 3rd round seeding [20]
    inputBinding:
      prefix: -y
  c:
    type: int?
    doc: skip seeds with more than INT occurrences [500]
    inputBinding:
      prefix: -c
  D:
    type: float?
    doc: |
      drop chains shorter than FLOAT fraction of the longest overlapping chain
      [0.50]
    inputBinding:
      prefix: -D
  W:
    type: int?
    doc: discard a chain if seeded bases shorter than INT [0]
    inputBinding:
      prefix: -W
  m:
    type: int?
    doc: perform at most INT rounds of mate rescues for each read [50]
    inputBinding:
      prefix: -m
  S:
    type: boolean
    default: false
    doc: skip mate rescue
    inputBinding:
      prefix: -S
  P:
    type: boolean
    default: false
    doc: skip pairing; mate rescue performed unless -S also in use
    inputBinding:
      prefix: -P

# Scoring options
  A:
    type: int?
    doc: |
      score for a sequence match, which scales options -TdBOELU unless
      overridden [1]
    inputBinding:
      prefix: -A
  B:
    type: int?
    doc: penalty for a mismatch [4]
    inputBinding:
      prefix: -B
  O:
    type: int[]?
    doc: gap open penalties for deletions and insertions [6,6]
    inputBinding:
      prefix: -O
      itemSeparator: ","
  E:
    type: int[]?
    doc: gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
    inputBinding:
      prefix: -E
      itemSeparator: ","
  L:
    type: int[]?
    doc: penalty for 5'- and 3'-end clipping [5,5]
    inputBinding:
      prefix: -L
      itemSeparator: ","
  U:
    type: int?
    doc: penalty for an unpaired read pair [17]
    inputBinding:
      prefix: -U
  x:
    type: string?
    doc: |
      read type. Setting -x changes multiple parameters unless overriden [null]
      pacbio -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 (PacBio reads to ref)
      ont2d -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0 (Oxford Nanopore 2D-reads to ref)
      intractg -B9 -O16 -L5 (intra-species contigs to ref)
    inputBinding:
      prefix: -x

# Input/output options
  p:
    type: boolean
    default: false
    doc: smart pairing (ignoring in2.fq)
    inputBinding:
      prefix: -p
  R:
    type: string?
    doc: read group header line such as '@RG\tID:foo\tSM:bar' [null]
    inputBinding:
      prefix: -R
  H:
    type: string?
    doc: |
      insert STR to header if it starts with @; or insert lines in FILE [null]
    inputBinding:
      prefix: -H
  j:
    type: boolean
    default: false
    doc: |
      treat ALT contigs as part of the primary assembly
      (i.e. ignore <idxbase>.alt file)
    inputBinding:
      prefix: -j
  v:
    type: int?
    doc: verbose level 1=error, 2=warning, 3=message, 4+=debugging [3]
    inputBinding:
      prefix: -v
  T:
    type: int?
    doc: minimum score to output [30]
    inputBinding:
      prefix: -T
  h:
    type: int[]?
    doc: |
      if there are <INT hits with score >80% of the max score, output all in XA
      [5,200]
    inputBinding:
      prefix: -h
  a:
    type: boolean
    default: false
    doc: output all alignments for SE or unpaired PE
    inputBinding:
      prefix: -a
  C:
    type: boolean
    default: false
    doc: append FASTA/FASTQ comment to SAM output
    inputBinding:
      prefix: -C
  V:
    type: boolean
    default: false
    doc: output the reference FASTA header in the XR tag
    inputBinding:
      prefix: -V
  Y:
    type: boolean
    default: false
    doc: use soft clipping for supplementary alignments
    inputBinding:
      prefix: -Y
  M:
    type: boolean
    doc: mark shorter split hits as secondary
    inputBinding:
      prefix: -M
  I:
    type: float[]?
    doc: |
      specify the mean, standard deviation (10% of the mean if absent), max
      (4 sigma from the mean if absent) and min of the insert size distribution.
      FR orientation only. [inferred]
    inputBinding:
      prefix: -I

outputs:
  output-file:
    type: File
    outputBinding:
      glob: $(inputs.output)

stdout: $(inputs.output)
