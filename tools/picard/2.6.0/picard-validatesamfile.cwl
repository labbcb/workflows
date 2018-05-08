cwlVersion: v1.0
class: CommandLineTool
id: picard-validatesamfile
label: "Picard ValidateSamFile"
doc: |
  Validates a SAM or BAM file. This tool reports on the validity of a SAM or BAM
  file relative to the SAM format specification. This is useful for
  troubleshooting errors encountered with other tools that may be caused by
  improper formatting, faulty alignments, incorrect flag values, etc. By
  default, the tool runs in VERBOSE mode and will exit after finding 100 errors
  and output them to the console (stdout). Therefore, it is often more practical
  to run this tool initially using the MODE=SUMMARY option.  This mode outputs a
  summary table listing the numbers of all 'errors' and 'warnings'.
  When fixing errors in your file, it is often useful to prioritize the severe
  validation errors and ignore the errors/warnings of lesser concern. This can
  be done using the IGNORE and/or IGNORE_WARNINGS arguments.
  After identifying and fixing your 'warnings/errors', we recommend that you
  rerun this tool to validate your SAM/BAM file prior to proceeding with your downstream analysis.  This will verify that all problems in your file have been addressed.

baseCommand: [java, -jar, /usr/picard.jar, ValidateSamFile]

hints:
  DockerRequirement:
    dockerPull: welliton/picard:2.6.0

inputs:
  input:
    type: File
    doc: Input SAM/BAM file.
    inputBinding:
      prefix: INPUT=
      separate: false
  output:
    type: string
    doc: Output file or standard out if missing.
    inputBinding:
      prefix: OUTPUT=
      separate: false
  mode:
    type: string?
    doc: |
      Mode of output  Default value: VERBOSE. This option can be set to 'null'
      to clear the default value. Possible values {VERBOSE, SUMMARY}
    inputBinding:
      prefix: MODE=
      separate: false
  ignore:
    type: string[]?
    doc: |
      List of validation error types to ignore.  Default value null. Possible
      values {INVALID_QUALITY_FORMAT, INVALID_FLAG_PROPER_PAIR,
      INVALID_FLAG_MATE_UNMAPPED, MISMATCH_FLAG_MATE_UNMAPPED,
      INVALID_FLAG_MATE_NEG_STRAND, MISMATCH_FLAG_MATE_NEG_STRAND,
      INVALID_FLAG_FIRST_OF_PAIR, INVALID_FLAG_SECOND_OF_PAIR,
      PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND,
      INVALID_FLAG_NOT_PRIM_ALIGNMENT, INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT,
      INVALID_FLAG_READ_UNMAPPED, INVALID_INSERT_SIZE, INVALID_MAPPING_QUALITY,
      INVALID_CIGAR, ADJACENT_INDEL_IN_CIGAR, INVALID_MATE_REF_INDEX,
      MISMATCH_MATE_REF_INDEX, INVALID_REFERENCE_INDEX, INVALID_ALIGNMENT_START,
      MISMATCH_MATE_ALIGNMENT_START, MATE_FIELD_MISMATCH, INVALID_TAG_NM,
      MISSING_TAG_NM, MISSING_HEADER, MISSING_SEQUENCE_DICTIONARY,
      MISSING_READ_GROUP, RECORD_OUT_OF_ORDER, READ_GROUP_NOT_FOUND,
      RECORD_MISSING_READ_GROUP, INVALID_INDEXING_BIN, MISSING_VERSION_NUMBER,
      INVALID_VERSION_NUMBER, TRUNCATED_FILE,
      MISMATCH_READ_LENGTH_AND_QUALS_LENGTH, EMPTY_READ,
      CIGAR_MAPS_OFF_REFERENCE, MISMATCH_READ_LENGTH_AND_E2_LENGTH,
      MISMATCH_READ_LENGTH_AND_U2_LENGTH, E2_BASE_EQUALS_PRIMARY_BASE,
      BAM_FILE_MISSING_TERMINATOR_BLOCK, UNRECOGNIZED_HEADER_TYPE,
      POORLY_FORMATTED_HEADER_TAG, HEADER_TAG_MULTIPLY_DEFINED,
      HEADER_RECORD_MISSING_REQUIRED_TAG, HEADER_TAG_NON_CONFORMING_VALUE,
      INVALID_DATE_STRING, TAG_VALUE_TOO_LARGE, INVALID_INDEX_FILE_POINTER,
      INVALID_PREDICTED_MEDIAN_INSERT_SIZE, DUPLICATE_READ_GROUP_ID,
      MISSING_PLATFORM_VALUE, INVALID_PLATFORM_VALUE,
      DUPLICATE_PROGRAM_GROUP_ID, MATE_NOT_FOUND, MATES_ARE_SAME_END,
      MISMATCH_MATE_CIGAR_STRING, MATE_CIGAR_STRING_INVALID_PRESENCE}
      This option may be specified 0 or more times.
    inputBinding:
      prefix: IGNORE=
      separate: false
  max_output:
    type: int?
    doc: |
      The maximum number of lines output in verbose mode Default value 100.
      This option can be set to 'null' to clear the default value.
  ignore_warnings:
    type: boolean
    default: false
    doc: |
      If true, only report errors and ignore warnings.  Default value: false.
      This option can be set to 'null' to clear the default value.
      Possible values {true, false}
    inputBinding:
      prefix: IGNORE_WARNINGS=true
  validate_index:
    type: string
    default: true
    doc: |
      DEPRECATED. Use INDEX_VALIDATION_STRINGENCY instead. If true and input is
      a BAM file with an index file, also validates the index.  Until this
      parameter is retired VALIDATE INDEX and INDEX_VALIDATION_STRINGENCY must
      agree on whether to validate the index. Default value: true. This option
      can be set to 'null' to clear the default value.
      Possible values {true, false}
    inputBinding:
      prefix: VALIDATE_INDEX=
      separate: false
  index_validation_stringency:
    type: int?
    doc: |
      If set to anything other than IndexValidationStringency. NONE and input
      is a BAM file with an index file, also validates the index at the
      specified stringency. Until VALIDATE_INDEX is retired, VALIDATE INDEX
      and INDEX_VALIDATION_STRINGENCY must agree on whether to validate the
      index. Default value EXHAUSTIVE. This option can be set to 'null' to
      clear the default value.
      Possible values {EXHAUSTIVE, LESS_EXHAUSTIVE, NONE}
    inputBinding:
      prefix: INDEX_VALIDATION_STRINGENCY=
      separate: false
  is_bisulfite_sequenced:
    type: boolean
    default: false
    doc: |
      Whether the SAM or BAM file consists of bisulfite sequenced reads. If so,
      C->T is not counted as an error in computing the value of the NM tag.
      Default value: false. This option can be set to 'null' to clear the
      default value. Possible values {true, false}
    inputBinding:
      prefix: IS_BISULFITE_SEQUENCED=true
  max_open_temp_files:
    type: int?
    doc: |
      Relevant for a coordinate-sorted file containing read pairs only. Maximum
      number of file handles to keep open when spilling mate info to disk. Set
      this number a little lower than the per-process maximum number of file
      that may be open. This number can be found by executing the 'ulimit -n'
      command on a Unix system.  Default value 8000. This option can be set to
      'null' to clear the default value.

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(inputs.output)
