task Faidx {
  File genomeFile

  command {
    ln -s ${genomeFile} -t .
    samtools faidx ${basename(genomeFile)}
  }

  output {
    File genomeIndexFile = "${basename(genomeFile)}.fai"
  }

  runtime {
    docker: "welliton/samtools:1.5"
  }
}
