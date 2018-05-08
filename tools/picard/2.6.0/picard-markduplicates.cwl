cwlVersion: v1.0
class: CommandLineTool
id: picard-markduplicates
label: "Picard MarkDuplicates"
doc: |
  Identifies duplicate reads. This tool locates and tags duplicate reads in a
  BAM or SAM file, where duplicate reads are defined as originating from a
  single fragment of DNA.

baseCommand: [java, -jar, /usr/picard.jar, MarkDuplicates]

hints:
  DockerRequirement:
    dockerPull: welliton/picard:2.6.0

inputs:
  input:
    type: File
    doc: |
      One or more input SAM or BAM files to analyze. Must be coordinate sorted.
    inputBinding:
      prefix: INPUT=
      separate: false
  output:
    type: string
    doc: The output file to write marked records to.
    inputBinding:
      prefix: OUTPUT=
      separate: false
  metrics_file:
    type: string?
    doc: File to write duplication metrics to.
    inputBinding:
      prefix: METRICS_FILE=
      separate: false
  remove_duplicates:
    type: boolean
    default: false
    doc: |
      If true do not write duplicates to the output file instead of writing them
      with appropriate flags set.  Default value false.
    inputBinding:
      prefix: REMOVE_DUPLICATES=true

  max_file_handles:
    type: int?
    doc: |
      Maximum number of file handles to keep open when spilling read ends to
      disk. Set this number a little lower than the per-process maximum number
      of file that may be open. This number can be found by executing the
      'ulimit -n' command on a Unix system.  Default value 8000.
    inputBinding:
      prefix: MAX_FILE_HANDLES=
      separate: false
  sorting_collection_size_ratio:
    type: float?
    doc: |
      This number, plus the maximum RAM available to the JVM, determine the
      memory footprint used by some of the sorting collections.  If you are
      running out of memory, try reducing this number.  Default value 0.25.
    inputBinding:
      prefix: SORTING_COLLECTION_SIZE_RATIO=
      separate: false
  barcode_tag:
    type: string?
    doc: Barcode SAM tag (ex. BC for 10X Genomics). Default value null.
    inputBinding:
      prefix: BARCODE_TAG=
      separate: false
  read_one_barcode:
    type: string?
    doc: |
      Read one barcode SAM tag (ex. BX for 10X Genomics) Default value null.
    inputBinding:
      prefix: READ_ONE_BARCODE_TAG=
      separate: false
  read_two_barcode:
    type: string?
    doc: Read two barcode SAM tag (ex. BX for 10X Genomics).
    inputBinding:
      prefix: READ_TWO_BARCODE_TAG
      separate: false
  tag_duplicate_set_members:
    type: boolean
    default: false
    doc: |
      If a read appears in a duplicate set, add two tags. The first tag,
      DUPLICATE_SET_SIZE_TAG (DS), indicates the size of the duplicate set.
      The smallest possible DS value is 2 which occurs when two reads map to
      the same portion of the reference only one of which is marked as
      duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a
      unique identifier for the duplicate set to which the record belongs.
      This identifier is the index-in-file of the representative read that was
      selected out of the duplicate set.  Default value false.
    inputBinding:
      prefix: TAG_DUPLICATE_SET_MEMBERS=true
  remove_sequencing_duplicates:
    type: boolean
    default: false
    doc: |
      If true remove 'optical' duplicates and other duplicates that appear to
      have arisen from the sequencing process instead of the library
      preparation process, even if REMOVE_DUPLICATES is false. If
      REMOVE_DUPLICATES is true, all duplicates are removed and this option is
      ignored. Default value false.
    inputBinding:
      prefix: REMOVE_SEQUENCING_DUPLICATES=true
  tagging_policy:
    type: string?
    doc: |
      Determines how duplicate types are recorded in the DT optional
      attribute.  Default value DontTag. This option can be set to 'null' to
      clear the default value. Possible values {DontTag, OpticalOnly, All}
    inputBinding:
      prefix: TAGGING_POLICY=
      separate: false
  assume_sort_order:
    type: string?
    doc: |
      If not null, assume that the input file has this order even if the
      header says otherwise. Possible values {unsorted, queryname, coordinate,
      duplicate, unknown}.
    inputBinding:
      prefix: ASSUME_SORT_ORDER=
      separate: false
  duplicate_scoring_strategy:
    type: string?
    doc: |
      The scoring strategy for choosing the non-duplicate among candidates.
      Default value SUM_OF_BASE_QUALITIES. This option can be set to 'null' to
      clear the default value. Possible values {SUM_OF_BASE_QUALITIES,
      TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}
    inputBinding:
      prefix: DUPLICATE_SCORING_STRATEGY=
      separate: false
  program_record_id:
    type: string?
    doc: |
      The program record ID for the @PG record(s) created by this program.
      Set to null to disable PG record creation.  This string may have a
      suffix appended to avoid collision with other program record IDs.
      Default value MarkDuplicates.
    inputBinding:
      prefix: PROGRAM_RECORD_ID=
      separate: false
  program_group_version:
    type: string?
    doc: |
      Value of VN tag of PG record to be created. If not specified, the
      version will be detected automatically.  Default value null.
    inputBinding:
      prefix: PROGRAM_GROUP_VERSION=
      separate: false
  program_group_command_line:
    type: string?
    doc: |
      Value of CL tag of PG record to be created. If not supplied the command
      line will be detected automatically.  Default value null.
    inputBinding:
      prefix: PROGRAM_GROUP_COMMAND_LINE=
      separate: false
  program_group_name:
    type: string?
    doc: |
      Value of PN tag of PG record to be created.
      Default value MarkDuplicates. This option can be set to 'null' to clear
      the default value.
    inputBinding:
      prefix: PROGRAM_GROUP_NAME=
      separate: false
  comment:
    type: string?
    doc: |
      Comment(s) to include in the output file's header. Default value null.
      This option may be specified 0 or more times.
    inputBinding:
      prefix: COMMENT=
      separate: false
  read_name_regex:
    type: string?
    doc: |
      Regular expression that can be used to parse read names in the incoming
      SAM file. Read names are parsed to extract three variables tile/region,
      x coordinate and y coordinate. These values are used to estimate the
      rate of optical duplication in order to give a more accurate estimated
      library size. Set this option to null to disable optical duplicate
      detection, e.g. for RNA-seq or other data where duplicate sets are
      extremely large and estimating library complexity is not an aim. Note
      that without optical duplicate counts, library size estimation will be
      inaccurate. The regular expression should contain three capture groups
      for the three variables, in order. It must match the entire read name.
      Note that if the default regex is specified, a regex match is not
      actually done, but instead the read name  is split on colon character.
      For 5 element names, the 3rd, 4th and 5th elements are assumed to be
      tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th,
      and 7th elements are assumed to be tile, x and y values.  Default value
       <optimized capture of last three ':' separated fields as numeric
       values>. This option can be set to 'null' to clear the default value.
    inputBinding:
      prefix: READ_NAME_REGEX=
      separate: false
  optical_duplicate_pixel_distance:
    type: string?
    doc: |
      The maximum offset between two duplicate clusters in order to consider
      them optical duplicates. The default is appropriate for unpatterned
      versions of the Illumina platform. For the patterned flowcell models,
      2500 is moreappropriate. For other platforms and models, users should
      experiment to find what works best.  Default value: 100. This option
      can be set to 'null' to clear the default value.
    inputBinding:
      prefix: OPTICAL_DUPLICATE_PIXEL_DISTANCE=
      separate: false
  max_optical_duplicate_set_size:
    type: int?
    doc: |
      This number is the maximum size of a set of duplicate reads for which we
      will attempt to determine which are optical duplicates.  Please be aware
      that if you raise this value too high and do encounter a very large set
      of duplicate reads, it will severely affect the runtime of this tool.
      To completely disable this check, set the value to -1.  Default value
      300000. This option can be set to 'null' to clear the default value.
    inputBinding:
      prefix: MAX_OPTICAL_DUPLICATE_SET_SIZE=
      separate: false

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output)
  output_metrics_file:
    type: File?
    outputBinding:
      glob: $(inputs.metrics_file)
