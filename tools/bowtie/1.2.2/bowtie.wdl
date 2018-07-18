task BuildIndex {
  Array[File] genomeFiles
  String ebwtBase
  
  command {
    bowtie-build ${sep=',' genomeFiles} ${ebwtBase}
  }
  
  output {
    Array[File] indexFiles = glob("*.ebwt")
  }
  
  runtime {
    docker: "welliton/bowtie:1.2.2"
  }
  
  parameter_meta {
    genomeFiles: "List of FASTA files containing the reference sequences to be aligned to."
    ebwtBase: "The basename of the index files to write."
  }
  
  meta {
    author: "Welliton Souza"
    email: "well309@gmail.com"
  }
}

task Align {
  Array[File] indexFiles
  String ebwtBase
  File file
  String outputFileName
  Int seedmms = 2
  Boolean all = false
  Boolean best = false
  Boolean strata = false
  Boolean sam = false
  Int threads = 1
  
  command {
    ln -s ${sep=' ' indexFiles} -t .
    bowtie \
      --seedmms ${seedmms} \
      ${true='--all' false='' all} \
      ${true='--best' false='' best} \
      ${true='--strata' false='' strata} \
      ${true='--sam' false='' strata} \
      --threads ${threads} \
      ${ebwtBase} \
      ${file} \
      ${outputFileName}
  }
  
  output {
    File outputFile = outputFileName
  }
  
  runtime {
    docker: "welliton/bowtie:1.2.2"
  }
  
  parameter_meta {
    indexFiles: "Index files (*.ebwt) created by BuidIndex task."
    ebwtBase: "The basename of the index to be searched."
    file: "A single file containing unpaired reads to be aligned."
    outputFileName: "File to write alignments to."
    seedmms: "Maximum number of mismatches permitted in the seed. This may be 0, 1, 2 or 3 and the default is 2."
    all: "Report all valid alignments per read or pair."
    best: "Make Bowtie guarantee that reported singleton alignments are best in terms of stratum."
    strata: "If many valid alignments exist and are reportable and they fall into more than one alignment stratum, report only those alignments that fall into the best stratum."
    sam: "Print alignments in SAM format."
    threads: "Launch parallel search threads (default: 1)."
  }
  
  meta {
    author: "Welliton Souza"
    email: "well309@gmail.com"
  }
}