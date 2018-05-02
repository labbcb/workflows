task GenomeGenerate {
  Array[File] genomeFiles
  File? gtfFile

  command {
    STAR --runMode --genomeGenerate --genomeDir . \
      --genomeFastaFiles ${sep=' ' genomeFiles} \
      ${'--sjdbGTFfile ' + gtfFile}
  }

  output {
    Array[File] indexFiles = glob("*")
  }

  runtime {
    docker: "welliton/star:2.5.3a"
  }
}

task AlignReads {
  Array[File] indexFiles
  Pair[File, File] pairedFiles

  command {
    mv ${sep=' ' indexFiles} .
    STAR --runMode alignReads --genomeDir . \
      --readFilesIn ${pairedFiles.left} ${pairedFiles.right} \
      --outSAMtype BAM SortedByCoordinate
  }

  output {
    File alignFile = "Aligned.sortedByCoord.out.bam"
    File logFinalFile = "Log.final.out"
    File logFile = "Log.out"
    File logProgressFile = "Log.progress.out"
    File sjFile = "SJ.out.tab"
  }

  runtime {
    docker: "welliton/star:2.5.3a"
  }
}
