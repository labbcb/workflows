task Index {
  File genomeFile
  String? algorithm

  command {
    ln -s ${genomeFile} -t .
    bwa index ${basename(genomeFile)} \
      ${'-a ' + algorithm}
  }

  output {
    File altFile = "${basename(genomeFile).alt"
    File saFile = "${basename(genomeFile).sa"
    File ambFile = "${basename(genomeFile).amb"
    File bwtFile = "${basename(genomeFile).bwt"
    File annFile = "${basename(genomeFile).ann"
    File pacFile = "${basename(genomeFile).pac"
  }

  runtime {
    docker: "welliton/bwa:0.7.15"
  }
}

task AlignMem {
  File genomeFile
  Array[File] indexFiles
  Pair[File, File] pairedFiles
  String outputFileName

  Boolean M = false
  String? R
  Int? t

  command {
    ln -s ${sep=' ' indexFiles} -t .
    bwa mem \
      ${true='-M' false='' M} \
      ${'-R ' + R} \
      ${'-t ' + t} \
      ${basename(genomeFile)} \
      ${pairedFiles.left} \
      ${pairedFiles.right} > ${outputFileName}
  }

  output {
    File alignFile = outputFileName
  }

  runtime {
    docker: "welliton/bwa:0.7.15"
  }
}
