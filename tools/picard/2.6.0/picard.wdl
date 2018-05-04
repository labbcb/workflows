task BuildBamIndex {
  File file

  command {
    mv ${file} .
    BuildBamIndex INPUT=${file}
  }

  output {
    File indexFile = sub(basename(file), ".bam$", ".bai")
  }

  runtime {
    docker: "welliton/picard:0.19.0"
  }
}

task CollectHsMetrics {
  File file
  String output
  File targetIntervals
  File baitIntervals

  File? referenceSequence
  File? referenceIndex
  String? perTargetCoverage

  command {
    CollectHsMetrics INPUT=${file} OUTPUT=${output} \
      TARGET_INTERVALS=${targetIntervals} \
      BAIT_INTERVALS=${baitIntervals} \
      REFERENCE_SEQUENCE=${referenceSequence} \
      ${'PER_TARGET_COVERAGE=' + perTargetCoverage}
  }

  output {
    File outputFile = output
    File? perTergetCoverageFile = perTargetCoverage
  }

  runtime {
    docker: "welliton/picard:0.19.0"
  }
}

task CreateSequenceDictionary {
  File reference
  String output

  String? genomeAssembly
  String? uri
  String? species
  String? truncateAtWhitespace
  Int? numSequences

  command {
    CreateSequenceDictionary \
      REFERENCE=${reference} \
      OUTPUT=${output} \
      ${'GENOME_ASSEMBLY=' + genomeAssembly} \

  }

  runtime {
    docker: "welliton/picard:0.19.0"
  }
}

task MarkDuplicates {
  File file
  String output
  String metricsFile?
  Boolean removeDuplicates = false

  command {
    MarkDuplicates INPUT=${file} \
      OUTPUT=${output} \
      ${'METRICS_FILE=' + metricsFile} \
      ${true='REMOVE_DUPLICATES=true' false='' removeDuplicates}
  }

  output {
    File outputFile = output
    File? metricsFile = metricsFile
  }

  runtime {
    docker: "welliton/picard:0.19.0"
  }
}

task SortSam {
  File file
  String output
  String? sortOrder

  command {
    SortSam INPUT=${file} \
      OUTPUT=${output} \
      ${'SORT_ORDER=' + sortOrder}
  }

  output {
    File outputFile = output
  }

  runtime {
    docker: "welliton/picard:0.19.0"
  }
}

task ValidateSamFile {
  File file
  String output
  String? mode

  command {
    ValidateSamFile INPUT=${input} \
      OUTPUT=${output} \
      ${'MODE=' + mode}
  }

  output {
    File outputFile = output
  }
}
