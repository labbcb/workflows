cwlVersion: v1.0
class: CommandLineTool
id: picard-sortsam
label: "Picard SortSam"
doc: |
  Sorts a SAM or BAM file.  This tool sorts the input SAM or BAM file by
  coordinate, queryname (QNAME), or some other property of the SAM record. The
  SortOrder of a SAM/BAM file is found in the SAM file header tag @HD in the
  field labeled SO.
  For a coordinate sorted SAM/BAM file, read alignments are sorted first by the
  reference sequence name (RNAME) field using the reference sequence dictionary
  (@SQ tag).  Alignments within these subgroups are secondarily sorted using the
  left-most mapping position of the read (POS).  Subsequent to this sorting
  scheme, alignments are listed arbitrarily.
  For queryname-sorted alignments, all alignments are grouped using the
  queryname field but the alignments are not necessarily sorted within these
  groups.  Reads having the same queryname are derived from the same template.

baseCommand: [java, -jar, /usr/picard.jar, SortSam]

hints:
  DockerRequirement:
    dockerPull: welliton/picard:2.6.0

inputs:
  input:
    type: File
    doc:  The BAM or SAM file to sort.
    inputBinding:
      prefix: INPUT=
      separate: false
  output:
    type: string
    doc: The sorted BAM or SAM output file.
    inputBinding:
      prefix: OUTPUT=
      separate: false
  sort-order:
    type: string?
    doc: |
      Sort order of output file.
      Possible values {unsorted, queryname, coordinate, duplicate, unknown}
    inputBinding:
      prefix: SORT_ORDER=
      separate: false

outputs:
  output-file:
    type: File
    outputBinding:
      glob: $(inputs.output)
