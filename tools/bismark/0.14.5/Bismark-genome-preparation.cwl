cwlVersion: v1.0
class: CommandLineTool

id: "bismark_genome_preparation"
label: "Bismark genome preparation"
doc: |
  Bismark is a program to map bisulfite treated sequencing reads to a genome of
  interest and perform methylation calls in a single step.

dct:creator:
  "@id": "https://orcid.org/0000-0001-8377-1836"
  foaf:name: Welliton Souza
  foaf:mbox: "mailto:well309@gmail.com"

baseCommand: [bismark_genome_preparation]

requirements:
  InitialWorkDirRequirement:
    listing: $(inputs.genome_files)
  DockerRequirement:
      dockerPull: welliton/bismark:0.14.5

arguments: ["."]

inputs:
  genome_files:
    type: File[]
    doc: |
      Genome files to be bisulfite converted. Bismark Genome Preparation expects
      one or more FASTA files (valid file extensions: '.fa' or '.fasta').
  path_to_bowtie:
    type: string?
    doc: |
      The full path to your Bowtie 1 or Bowtie 2 installation. Only required if
      Bowtie/Bowtie2 is not in the PATH.
    inputBinding:
      prefix: --path_to_bowtie
  bowtie1:
    type: boolean?
    doc: This will create bisulfite indexes for Bowtie 1. Default Bowtie 2.
    inputBinding:
      prefix: --bowtie1
  bowtie2:
    type: boolean?
    doc: This will create bisulfite indexes for Bowtie 2. Default ON.
    inputBinding:
      prefix: --bowtie2
  single_fasta:
    type: boolean
    doc: |
      Instruct the Bismark Indexer to write the converted genomes into
      single-entry FastA files instead of making one multi-FastA file (MFA) per
      chromosome. This might be useful if individual bisulfite converted
      chromosomes are needed (e.g. for debugging), however it can cause a
      problem with indexing if the number of chromosomes is vast (this is likely
      to be in the range of several thousand files; operating systems can only
      handle lists up to a certain length. Some newly assembled genomes may
      contain 20000-50000 contig of scaffold files which do exceed this list
      length limit).
    default: false
    inputBinding:
      prefix: --single_fasta
  genomic_composition:
    type: boolean
    doc: |
      Calculate and extract the genomic sequence composition for mono- and
      di-nucleotides and write the genomic composition table
      genomic_nucleotide_frequencies.txt to the genome folder. This may be
      useful later on when using bam2nuc or the Bismark option
      'nucleotide_coverage'.
    default: false
    inputBinding:
      prefix: --genomic_composition

outputs:
  genome_folder:
    type: Directory
    outputBinding:
      glob: "Bisulfite_Genome"
