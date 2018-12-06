version 1.0

task Faidx {

    input {
      File genomeFile
    }

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
