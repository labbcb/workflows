task Faidx {
  File genomeFile

  command {
    samtools faidx ${genomeFile}
  }

  output {
    File genomeIndexFile = "${basename(genomeFile)}.fai"
  }

  runtime {
    docker: "welliton/samtools:1.5"
  }
}
