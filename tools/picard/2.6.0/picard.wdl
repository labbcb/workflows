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
    docker: "welliton/picard:2.6.0"
  }
}

task CollectHsMetrics {
  File file
  File indexFile
  String outputFileName
  File targetIntervalsFile
  File baitIntervalsFile
  String verbosity = "INFO"

  File? genomeFile
  File? referenceIndex
  String? perTargetCoverage

  command {
    CollectHsMetrics INPUT=${file} \
      OUTPUT=${outputFileName} \
      TARGET_INTERVALS=${targetIntervalsFile} \
      BAIT_INTERVALS=${baitIntervalsFile} \
      REFERENCE_SEQUENCE=${genomeFile} \
      ${'PER_TARGET_COVERAGE=' + perTargetCoverage} \
      ${true='VERBOSITY=true' false='' verbosity}
  }

  output {
    File outputFile = outputFileName
    File? perTergetCoverageFile = perTargetCoverage
  }

  runtime {
    docker: "welliton/picard:2.6.0"
  }
}

task CreateSequenceDictionary {
  File genomeFile
  String outputFileName

  String? genomeAssembly
  String? uri
  String? species
  Boolean? truncateNamesAtWhitespace
  Int? numSequences

  command {
    CreateSequenceDictionary \
      REFERENCE=${genomeFile} \
      OUTPUT=${outputFileName} \
      ${'GENOME_ASSEMBLY=' + genomeAssembly} \
      ${'URI=' + uri} \
      ${'SPECIES=' + species} \
      ${true='TRUNCATE_NAMES_AT_WHITESPACE=true' false='' truncateNamesAtWhitespace} \
      ${'NUM_SEQUENCES=' + numSequences}
  }

  output {
    File dictionaryFile = outputFileName
  }

  runtime {
    docker: "welliton/picard:2.6.0"
  }
}

task MarkDuplicates {
  File file
  String outputFileName
  String? metricsFileName
  Boolean removeDuplicates = false

  command {
    MarkDuplicates INPUT=${file} \
      OUTPUT=${outputFileName} \
      ${'METRICS_FILE=' + metricsFileName} \
      ${true='REMOVE_DUPLICATES=true' false='' removeDuplicates}
  }

  output {
    File outputFile = outputFileName
    File? metricsFile = metricsFileName
  }

  runtime {
    docker: "welliton/picard:2.6.0"
  }
}

task SortSam {
  File file
  String outputFileName
  String? sortOrder

  command {
    SortSam INPUT=${file} \
      OUTPUT=${outputFileName} \
      ${'SORT_ORDER=' + sortOrder}
  }

  output {
    File outputFile = outputFileName
  }

  runtime {
    docker: "welliton/picard:2.6.0"
  }
}

task ValidateSamFile {
  File file
  File indexFile
  String outputFileName
  String? mode
  Boolean validateIndex = true

  command {
    ValidateSamFile INPUT=${file} \
      OUTPUT=${outputFileName} \
      ${'MODE=' + mode} \
      VALIDATE_INDEX=${validateIndex}
  }

  output {
    File outputFile = outputFileName
  }
}
