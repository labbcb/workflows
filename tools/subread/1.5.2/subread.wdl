version 1.0

task FeatureCounts {

    input {
      File annotationFile
      Array[File] files
      String outputFileName

      String format = "GTF"
      String featureType = "exon"
      String attributeType = "gene_id"
      Boolean featureLevel = false
      Boolean overlaps = false
      Boolean multiMapping = false
      Int mappingQuality = 0
      Int threads = 1
    }
  
  command {
    featureCounts \
      -F ${format} \
      -t ${featureType} \
      -g ${attributeType} \
      ${true='-f' false='' featureLevel} \
      ${true='-O' false='' overlaps} \
      ${true='-M' false='' multiMapping} \
      -Q ${mappingQuality} \
      -T ${threads} \
      -a ${annotationFile} \
      -o ${outputFileName} \
      ${sep=' ' files}
  }
  
  output {
    File outputFile = outputFileName
    File summaryFile = "${outputFileName}.summary"
  }

  runtime {
    docker: "welliton/subread:1.5.2"
  }
}