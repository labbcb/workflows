task BuildBamIndex {
  File file

  command {
    ln -s ${file} -t .
    java -jar /usr/picard.jar BuildBamIndex INPUT=${basename(file)}
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

  File genomeFile
  File genomeIndexFile
  String? perTargetCoverage

  command {
    ln -s ${file} ${indexFile} ${genomeFile} ${genomeIndexFile} -t .
    java -jar /usr/picard.jar CollectHsMetrics \
      INPUT=${basename(file)} \
      OUTPUT=${outputFileName} \
      TARGET_INTERVALS=${targetIntervalsFile} \
      BAIT_INTERVALS=${baitIntervalsFile} \
      REFERENCE_SEQUENCE=${basename(genomeFile)} \
      ${'PER_TARGET_COVERAGE=' + perTargetCoverage} \
      ${'VERBOSITY=' + verbosity}
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
    java -jar /usr/picard.jar CreateSequenceDictionary \
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
    java -jar /usr/picard.jar MarkDuplicates INPUT=${file} \
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
    java -jar /usr/picard.jar SortSam INPUT=${file} \
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
    ln -s ${file} ${indexFile} -t .
    java -jar /usr/picard.jar ValidateSamFile \
      INPUT=${basename(file)} \
      OUTPUT=${outputFileName} \
      ${'MODE=' + mode} \
      VALIDATE_INDEX=${validateIndex}
  }

  output {
    File outputFile = outputFileName
  }

  runtime {
    docker: "welliton/picard:2.6.0"
  }
}

task FastqToSam {
  File fileR1
  File? fileR2
  Boolean useSequentialFastqs = false
  String? qualityFormat
  String outputFileName
  String? readGroupName
  String sampleName
  String? libraryName
  String? platformUnit
  String? platform
  String? sequencingCenter
  Int? predictedInsertSize
  String? programGroup
  String? platformModel
  String? comment
  String? description
  String? runDate
  String sortOrder = "queryname"
  Int minQ = 0
  Int maxQ = 93
  Boolean allowAndIgnoreEmptyLines = false

  command {
    java -jar /usr/picard.jar FastqToSam \
      FASTQ=${fileR1} \
      ${'FASTQ2=' + fileR2} \
      ${'USE_SEQUENTIAL_FASTQS=' + useSequentialFastqs} \
      ${'QUALITY_FORMAT=' + qualityFormat} \
      OUTPUT=${outputFileName} \
      ${'READ_GROUP_NAME=' + readGroupName} \
      SAMPLE_NAME=${sampleName} \
      ${'LIBRARY_NAME=' + libraryName} \
      ${'PLATFORM_UNIT=' + platformUnit} \
      ${'PLATFORM=' + platform} \
      ${'SEQUENCING_CENTER=' + sequencingCenter} \
      ${'PREDICTED_INSERT_SIZE=' + predictedInsertSize} \
      ${'PROGRAM_GROUP=' + programGroup} \
      ${'PLATFORM_MODEL=' + platformModel} \
      ${'COMMENT=' + comment} \
      ${'DESCRIPTION=' + description} \
      ${'RUN_DATE=' + runDate} \
      SORT_ORDER=${sortOrder} \
      MIN_Q=${minQ} \
      MAX_Q=${maxQ} \
      ALLOW_AND_IGNORE_EMPTY_LINES=${allowAndIgnoreEmptyLines}
  }

  output {
    File outputFile = outputFileName
  }

  runtime {
    docker: "welliton/picard:2.6.0"
  }

}
