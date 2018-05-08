task Index {
  File genomeFile
  String? algorithm

  command {
    bwa index ${genomeFile} \
      ${'-a ' + algorithm}
  }

  output {
    Array[File] indexFiles = glob("*")
  }

  runtime {
    docker: "welliton/bwa:0.7.15"
  }
}

task AlignMem {
  File genomeFile
  Array[File] indexFiles
  Pair[File, File] pairedFiles

  Boolean M = false
  String? R
  Int? t

  command {
    mv ${sep=' ' indexFiles} .
    bwa mem ${basename(genomeFile)} \
      ${pairedFiles.left} ${pairedFiles.right} \
      ${true='-M' false='' M} \
      ${'-R ' + R} \
      ${'-t ' + t}
  }

  output {
    File alignFile = stdout()
  }
}
