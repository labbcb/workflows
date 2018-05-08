task Faidx {
  File genomeFile

  command {
    samtools faidx ${genomeFile}
  }

  output {
    File genomeIndexFile = "${genomeFile}.fai"
  }

  runtime {
    docker: "welliton/samtools:1.5"
  }
}
