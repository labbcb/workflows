task Index {
  File genomeFile
  String? algorithm

  command {
    ln -s ${genomeFile} -t .
    bwa index ${basename(genomeFile)} \
      ${'-a ' + algorithm}
  }

  output {
    Array[File] indexFiles = glob("${basename(genomeFile)}.*")
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
    bwa mem ${basename(genomeFile)} \
      ${pairedFiles.left} ${pairedFiles.right} \
      ${true='-M' false='' M} \
      ${'-R ' + R} \
      ${'-t ' + t} > ${outputFileName}
  }

  output {
    File alignFile = outputFileName
  }

  runtime {
    docker: "welliton/bwa:0.7.15"
  }
}
