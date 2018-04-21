cwlVersion: v1.0
class: CommandLineTool

id: "Bismark"
label: "bismark"
doc: |
  Bismark is a program to map bisulfite treated sequencing reads to a genome of
  interest and perform methylation calls in a single step.

dct:creator:
  "@id": "https://orcid.org/0000-0001-8377-1836"
  foaf:name: Welliton Souza
  foaf:mbox: "mailto:well309@gmail.com"

baseCommand: [bismark]

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: |
      ${
        var listing = inputs.genome_files;
        listing.push(inputs.genome_folder);
        return listing;
       }
  DockerRequirement:
      dockerPull: welliton/bismark:v0.14.5

arguments: ["."]

inputs:

# Mandatory
  genome_folder:
    type: Directory
    doc: |
      Directory 'Bisulfite_Genome' created by 'bismark_genome_preparation'.
  genome_files:
    type: File[]
    doc: Unmodified reference genome files used by 'bismark_genome_preparation'.

  mate1:
    type: File?
    doc: A single file containing the mate number 1.
    inputBinding:
      prefix: "-1"
      position: 2
  mate2:
    type: File?
    doc: A single file containing the mate number 2.
    inputBinding:
      prefix: "-2"
      position: 2
  single:
    type: File?
    doc: |
      A single file containing the reads to be aligned.
    inputBinding: {}

# Input
  single_end:
    type: boolean
    doc: Sets single-end mapping mode explicitly giving a list of file names.
    default: false
    inputBinding:
      prefix: --single_end
  fastq:
    type: boolean
    doc: The query input files are FASTQ files. Default.
    default: true
    inputBinding:
      prefix: --fastq
  fasta:
    type: boolean
    doc: |
      The query input files are FASTA files. All quality values are assumed to
      be 40 on the Phred scale. FASTA files are expected to contain the read
      name and the sequence on a single line each (and not spread over several
      lines).
    default: false
    inputBinding:
      prefix: --fasta
  skip:
    type: int?
    doc: |
      Skip (i.e. do not align) the first number of reads or read pairs from the
      input.
    inputBinding:
      prefix: --skip
  upto:
    type: int?
    doc: |
      Only aligns the first number of reads or read pairs from the input.
      Default no limit.
    inputBinding:
      prefix: --upto
  phred33_quals:
    type: boolean
    doc: |
      FastQ qualities are ASCII chars equal to the Phred quality plus 33.
      Default ON.
    default: true
    inputBinding:
      prefix: --phred33-quals
  phred64_quals:
    type: boolean
    doc: |
      FastQ qualities are ASCII chars equal to the Phred quality plus 64.
      Default OFF.
    default: false
    inputBinding:
      prefix: --phred64-quals
  solexa_quals:
    type: boolean
    doc: |
      Convert FastQ qualities from solexa-scaled (which can be negative) to
      phred-scaled (which can't). The formula for conversion is:
      phred-qual = 10 * log(1 + 10 ** (solexa-qual/10.0)) / log(10).
      Used with -q. This is usually the right option for use with (unconverted)
      reads emitted by the GA Pipeline versions prior to 1.3. Default OFF.
    default: false
    inputBinding:
      prefix: --solexa1.3-quals
  solexa1_3_quals:
    type: boolean
    doc: |
      Same as phred64_quals. This is usually the right option for use with
      (unconverted) reads emitted by GA Pipeline version 1.3 or later.
      Default OFF.
    default: false
    inputBinding:
      prefix: --solexa1.3-quals

# Alignment

  seedmms:
    type: int
    doc: |
      The maximum number of mismatches permitted in the "seed", i.e. the first L
      base pairs of the read (where L is set with seedlen). This may be 0, 1, 2
      or 3 and the default is 1. This option is only available for Bowtie 1.
    default: 1
    inputBinding:
      prefix: --seedmms
  seedlen:
    type: int
    doc: |
      The "seed length"; i.e., the number of bases of the high quality end of
      the read to which the -n ceiling applies. The default is 28. Bowtie (and
      thus Bismark) is faster for larger values of -l. This option is only
      available for Bowtie 1.
    default: 28
    inputBinding:
      prefix: --seedlen
  maqerr:
    type: int
    doc: |
      Maximum permitted total of quality values at all mismatched read positions
      throughout the entire alignment, not just in the "seed".
      The default is 70. Like Maq, Bowtie rounds quality values to the nearest
      10 and saturates at 30.
    default: 70
    inputBinding:
      prefix: --maqerr
  chunkmbs:
    type: int
    doc: |
      The number of megabytes of memory a given thread is given to store path
      descriptors in --best mode. Best-first search must keep track of many
      paths at once to ensure it is always extending the path with the lowest
      cumulative cost. Bowtie tries to minimize the memory impact of the
      descriptors, but they can still grow very large in some cases. If you
      receive an error message saying that chunk memory has been exhausted in
      best mode, try adjusting this parameter up to dedicate more memory to
      the descriptors. Default 512.
    default: 512
    inputBinding:
      prefix: --chunkmbs
  minins:
    type: int
    doc: |
      The minimum insert size for valid paired-end alignments. E.g. if -I 60 is
      specified and a paired-end alignment consists of two 20-bp alignments in
      the appropriate orientation with a 20-bp gap between them, that alignment
      is considered valid (as long as maxins is also satisfied).
      A 19-bp gap would not be valid in that case. Default 0.
    default: 0
  maxins:
    type: int
    doc: |
      The maximum insert size for valid paired-end alignments. E.g. if -X 100 is
      specified and a paired-end alignment consists of two 20-bp alignments in
      the proper orientation with a 60-bp gap between them, that alignment is
      considered valid (as long as minins is also satisfied).
      A 61-bp gap would not be valid in that case. Default 500.
    default: 500
  multicore:
    type: int?
    doc: |
      Sets the number of parallel instances of Bismark to be run concurrently.
      This forks the Bismark alignment step very early on so that each individu
      al Spawn of Bismark processes only every n-th sequence (N being set by
      multicore). Once all processes have completed, the individual BAM files,
      mapping reports, unmapped or ambiguous FastQ files are merged into single
      files in very much the same way as they would have been generated running
      Bismark conventionally with only a single instance.

# Bowtie 1 Reporting

  k:
    type: int
    doc: |
      Due to the way Bismark works Bowtie 1 will report up to 2 valid
      alignments. This option is used by default and cannot be changed.
    default: 2
  no_best:
    type: boolean
    doc: |
      Disables the 'best' option which is on by default. This can speed up the
      alignment process, e.g. for testing purposes, but for credible results it
      is not recommended to disable 'best'.
    default: false
    inputBinding:
      prefix: --no_best

# Output

  non_directional:
    type: boolean
    doc: |
      The sequencing library was constructed in a non strand-specific manner,
      alignments to all four bisulfite strands will be reported. (The current
      Illumina protocol for BS-Seq is directional, in which case the strands
      complementary to the original strands are merely theoretical and should
      not exist in reality. Specifying directional alignments (which is the
      default) will only run 2 alignment threads to the original top (OT) or
      bottom (OB) strands in parallel and report these alignments. This is the
      recommended option for strand-specific libraries). Default OFF.
    default: false
    inputBinding:
      prefix: --non_directional
  pbat:
    type: boolean
    doc: |
      This option may be used for PBAT-Seq libraries (Post-Bisulfite Adapter
      Tagging; Kobayashi et al., PLoS Genetics, 2012). This is essentially the
      exact opposite of alignments in 'directional' mode, as it will only launch
      two alignment threads to the CTOT and CTOB strands instead of the normal
      OT and OB ones. Use this option only if you are certain that your
      libraries were constructed following a PBAT protocol (if you don't know
      what PBAT-Seq is you should not specify this option). The option 'pbat'
      works only for FastQ files and uncompressed temporary files.
    default: false
    inputBinding:
      prefix: --pbat
  sam_no_hd:
    type: boolean
    doc: |
      Suppress SAM header lines (starting with @). This might be useful when
      very large input files are split up into several smaller files to run
      concurrently and the output files are to be merged afterwards.
    default: false
    inputBinding:
      prefix: --sam-no-hd
  rg_tag:
    type: boolean
    doc: |
      Write out a Read Group tag to the resulting SAM/BAM file. This will write
      the following line to the SAM header:
        @RG PL: ILLUMINA ID:SAMPLE SM:SAMPLE
      to set ID and SM see 'rg_id' and 'rg_sample'. In addition each read
      receives an RG:Z:RG-ID tag. Default OFF (to not inflate file sizes).
    default: false
    inputBinding:
      prefix: --rg_tag
  rg_id:
    type: string?
    doc: Sets the ID field in the @RG header line. Default SAMPLE
    inputBinding:
      prefix: --rg_id
  rg_sample:
    type: string?
    doc: |
      Sets the SM field in the @RG header line; can't be set without setting
      'rg_id' as well. Default SAMPLE
    inputBinding:
      prefix: --rg_sample
  quiet:
    type: boolean
    doc: Print nothing besides alignments.
    default: false
    inputBinding:
      prefix: --quiet
  vanilla:
    type: boolean
    doc: |
      Performs bisulfite mapping with Bowtie 1 and prints the 'old' custom
      Bismark output (up to versions 0.5.X) instead of SAM format output.
    default: false
    inputBinding:
      prefix: --vanilla
  un:
    type: boolean
    doc: |
      Write all reads that could not be aligned to the file
      '_unmapped_reads.fq.gz' in the output directory. Written reads will appear
      as they did in the input, without any translation of quality values that
      may have taken place within Bowtie or Bismark. Paired-end reads will be
      written to two parallel files with _1 and _2 inserted in their filenames,
      i.e. unmapped_reads_1.fq.gz and unmapped_reads_2.fq.gz. Reads with more
      than one valid alignment with the same number of lowest mismatches
      (ambiguous mapping) are also written to 'unmapped_reads.fq.gz' unless
      'ambiguous' is also specified.
    default: false
    inputBinding:
      prefix: --un
  ambiguous:
    type: boolean
    doc: |
      Write all reads which produce more than one valid alignment with the same
      number of lowest mismatches or other reads that fail to align uniquely to
      '_ambiguous_reads.fq'. Written reads will appear as they did in the input,
      without any of the translation of quality values that may have taken place
      within Bowtie or Bismark. Paired-end reads will be written to two parallel
      files with '_1' and '_2' inserted in their filenames, i.e.
      '_ambiguous_reads_1.fq' and '_ambiguous_reads_2.fq'. These reads are not
      written to the file specified with 'un'.
    default: false
    inputBinding:
      prefix: --ambigous
  output_dir:
    type: string?
    doc: |
      Write all output files into this directory. By default the output files
      will be written into the same folder as the input file. If the specified
      folder does not exist, Bismark will attempt to create it first. The path
      to the output folder can be either relative or absolute.
    inputBinding:
      prefix: --output_dir
  temp_dir:
    type: string?
    doc: |
      Write temporary files to this directory instead of into the same directory
      as the input files. If the specified folder does not exist, Bismark will
      attempt to create it first. The path to the temporary folder can be either
      relative or absolute.
    inputBinding:
      prefix: --temp_dir
  non_bs_mm:
    type: boolean
    doc: |
      Optionally outputs an extra column specifying the number of non-bisulfite
      mismatches a read during the alignment step. This option is only available
      for SAM format. In Bowtie 2 context, this value is just the number of
      actual non-bisulfite mismatches and ignores potential insertions or
      deletions. The format for single-end reads and read 1 of paired-end reads
      is 'XA:Z:number of mismatches' and 'XB:Z:number of mismatches' for read 2
      of paired-end reads.
    default: false
    inputBinding:
      prefix: --non_bs_mm
  gzip:
    type: boolean
    doc: |
      Temporary bisulfite conversion files will be written out in a GZIP
      compressed form to save disk space. This option is available for most
      alignment modes but is not available for paired-end FASTA files.
    default: false
    inputBinding:
      prefix: --gzip
  sam:
    type: boolean
    doc: |
      The output will be written out in SAM format instead of the default BAM
      format.
    default: false
    inputBinding:
      prefix: --sam
  cram:
    type: boolean
    doc: |
      Writes the output to a CRAM file instead of BAM. This requires the use of
      Samtools 1.2 or higher.
    default: false
    inputBinding:
      prefix: --cram
  cram_ref:
    type: File?
    doc: |
      CRAM output requires you to specify a reference genome as a single FASTA
      file. If this single-FastA reference file is not supplied explicitly it
      will be regenerated from the genome .fa sequence(s) used for the Bismark
      run and written to a file called 'Bismark_genome_CRAM_reference.mfa' into
      the output directory.
    inputBinding:
      prefix: --cram_ref
  samtools_path:
    type: string?
    doc: |
      The path to your Samtools installation, e.g. '/home/user/samtools/'. Does
      not need to be specified explicitly if Samtools is in the PATH already.
    inputBinding:
      prefix: --samtools_path
  prefix:
    type: string?
    doc: |
      Prefixes 'prefix' to the output file names. Trailing dots will be
      replaced by a single one. For example, 'prefix test' with 'file.fq'
      would result in the output file 'test.file_bismark.bam' etc.
    inputBinding:
      prefix: --prefix
  basename:
    type: string?
    doc: |
      Write all output to files starting with this base file name. For example,
      'basename foo' would result in the files 'foo.bam' and 'foo_SE_report.txt'
      (or its paired-end equivalent). Takes precedence over 'prefix'.
    inputBinding:
      prefix: --basename
  old_flag:
    type: boolean
    doc: |
      Only in paired-end SAM mode, uses the FLAG values used by Bismark 0.8.2
      and before. In addition, this options appends /1 and /2 to the read IDs
      for reads 1 and 2 relative to the input file. Since both the appended read
      IDs and custom FLAG values may cause problems with some downstream tools
      such as Picard, new defaults were implemented as of version 0.8.3.
    default: false
    inputBinding:
      prefix: --old_flag
  ambig_bam:
    type: boolean
    doc: |
      For reads that have multiple alignments a random alignment is written out
      to a special file ending in '.ambiguous.bam'. The alignments are in
      Bowtie2 format and do not any contain Bismark specific entries such as the
      methylation call etc. These ambiguous BAM files are intended to be used as
      coverage estimators for variant callers.
    default: false
    inputBinding:
      prefix: --ambig_bam
  nucleotide_coverage:
    type: boolean
    doc: |
      Calculates the mono- and di-nucleotide sequence composition of covered
      positions in the analysed BAM file and compares it to the genomic average
      composition once alignments are complete by calling  bam2nuc. Since this
      calculation may take a while, bam2nuc attempts to write the genomic
      sequence composition into a file called
      'genomic_nucleotide_frequencies.txt' inside the reference genome folder so
      it can be re-used the next time round instead of calculating it once
      again. If a file nucleotide_stats.txt is found with the Bismark reports it
      will be automatically detected and used for the Bismark HTML report. This
      option works only for BAM or CRAM files.
    default: false
    inputBinding:
      prefix: --nucleotide_coverage

# Other

  bowtie1:
    type: boolean
    doc: |
      Uses Bowtie 1 instead of Bowtie 2, which might be a good choice for faster
      and very short alignments. Bismark assumes that raw sequence data is
      adapter and/or quality trimmed where appropriate. Both small (.ebwt) and
      large (.ebwtl) Bowtie indexes are supported. Default OFF.
    default: false
    inputBinding:
      prefix: --bowtie1

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bam"
  report:
    type: File
    outputBinding:
      glob: "*_report.txt"
