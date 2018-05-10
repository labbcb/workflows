task BaseRecalibrator {
  File file
  File indexFile
  File genomeFile
  File genomeIndexFile
  File dictionaryFile
  String outputFileName
  File? intervalsFile
  Int? intervalPadding
  Array[File] knownSitesFiles
  Array[File] knownSitesIndexFiles
  File? BQSR

  command {
    ln -s ${file} ${indexFile} ${genomeFile} ${genomeIndexFile} ${dictionaryFile} -t .
    java -jar /usr/GenomeAnalysisTK.jar --analysis_type BaseRecalibrator \
      --input_file ${basename(file)} \
      --reference_sequence ${basename(genomeFile)} \
      --out ${outputFileName} \
      --knownSites ${sep=' --knownSites ' knownSitesFiles} \
      ${'--intervals ' + intervalsFile} \
      ${'--interval_padding ' + intervalPadding} \
      ${'--BQSR ' + BQSR}
  }

  output {
    File outputFile = outputFileName
  }

  runtime {
    docker: "broadinstitute/gatk3:3.6-0"
  }
}

task AnalyzeCovariates {
  File beforeReportFile
  File afterReportFile
  File genomeFile
  File genomeIndexFile
  File dictionaryFile
  String outputFileName
  File? intervalsFile
  Int? intervalPadding

  command {
    ln -s ${genomeFile} ${genomeIndexFile} ${dictionaryFile} -t .
    java -jar /usr/GenomeAnalysisTK.jar --analysis_type AnalyzeCovariates \
      --beforeReportFile ${beforeReportFile} \
      --afterReportFile ${afterReportFile} \
      --reference_sequence ${basename(genomeFile)} \
      --plotsReportFile ${outputFileName} \
      ${'--intervals ' + intervalsFile} \
      ${'--interval_padding ' + intervalPadding}
  }

  output {
    File outputFile = outputFileName
  }

  runtime {
    docker: "broadinstitute/gatk3:3.6-0"
  }
}

task PrintReads {
  File file
  File indexFile
  File genomeFile
  File genomeIndexFile
  File dictionaryFile
  String outputFileName
  File? intervalsFile
  Int? intervalPadding
  File? BQSR

  command {
    ln -s ${file} ${indexFile} ${genomeFile} ${genomeIndexFile} ${dictionaryFile} -t .
    java -jar /usr/GenomeAnalysisTK.jar --analysis_type PrintReads \
      --input_file ${basename(file)} \
      --reference_sequence ${basename(genomeFile)} \
      --out ${outputFileName} \
      ${'--intervals ' + intervalsFile} \
      ${'--interval_padding ' + intervalPadding} \
      ${'--BQSR ' + BQSR}
  }

  output {
    File outputFile = outputFileName
    File outputIndexFile = sub(outputFileName, ".bam$", ".bai")
  }

  runtime {
    docker: "broadinstitute/gatk3:3.6-0"
  }
}

task HaplotypeCaller {
  File file
  File indexFile
  File genomeFile
  File genomeIndexFile
  File dictionaryFile
  String outputFileName
  File? intervalsFile
  Int? intervalPadding
  String? genotypingMode
  Int? standCallConf
  Int? standEmitConf
  String? emitRefConfidence
  String? variantIndexType
  Int? variantIndexParameter

  command {
    ln -s ${file} ${indexFile} ${genomeFile} ${genomeIndexFile} ${dictionaryFile} -t .
    java -jar /usr/GenomeAnalysisTK.jar --analysis_type HaplotypeCaller \
      --input_file ${file} \
      --reference_sequence ${basename(genomeFile)} \
      --out ${outputFileName} \
      ${'--intervals ' + intervalsFile} \
      ${'--interval_padding ' + intervalPadding} \
      ${'--genotyping_mode ' + genotypingMode} \
      ${'-stand_call_conf ' + standCallConf} \
      ${'-stand_emit_conf ' + standEmitConf} \
      ${'--emitRefConfidence ' + emitRefConfidence} \
      ${'--variant_index_type ' + variantIndexType} \
      ${'--variant_index_parameter ' + variantIndexParameter}
  }

  output {
    File vcfFile = outputFileName
    File vcfIndexFile = "${outputFileName}.idx"
  }

  runtime {
    docker: "broadinstitute/gatk3:3.6-0"
  }
}
