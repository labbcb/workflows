import "https://raw.githubusercontent.com/labbcb/rnnr/master/tools/trimgalore/0.5.0/trimgalore.wdl" as trimgalore
import "https://raw.githubusercontent.com/labbcb/rnnr/master/tools/bowtie/1.2.2/bowtie.wdl" as bowtie
import "https://raw.githubusercontent.com/labbcb/rnnr/master/tools/subread/1.6.2/subread.wdl" as subread

workflow MirnaSeq {
  
  Array[File] files
  Array[File] genomeFiles
  String ebwtBase
  File annotationFile
  
  String countFileName = "counts.txt"
  Int workers = 1
  Boolean smallRNA = true
  Int length = 18
  Int maxN = 2
  Int maxLength = 30
  Boolean gzip = false
  Int seedmms = 2
  Boolean all = true
  Boolean best = true
  Boolean strata = true
  Boolean sam = true
  String featureType = "miRNA"
  String attributeType = "ID"
  Boolean overlaps = true
  Boolean multiMapping = true
  Int mappingQuality = 20
  
  call bowtie.BuildIndex {
    input:
      genomeFiles = genomeFiles,
      ebwtBase = ebwtBase
  }
  
  scatter (file in files) {
    call trimgalore.TrimGaloreSingle {
      input:
        file = file,
        smallRNA = smallRNA,
        length =  length,
        maxLength = maxLength,
        maxN = maxN,
        gzip = gzip
    }
    call bowtie.Align {
      input:
        indexFiles = BuildIndex.indexFiles,
        ebwtBase = ebwtBase,
        file = TrimGaloreSingle.trimFile,
        outputFileName = sub(basename(file), "_R1.fastq.gz", ".sam"),
        seedmms = seedmms,
        all = all,
        best = best,
        strata = strata,
        sam = sam,
        threads = workers
    }
  }
  
  call subread.FeatureCounts {
    input:
      files = Align.outputFile,
      annotationFile = annotationFile,
      featureType = featureType,
      attributeType = attributeType,
      overlaps = overlaps,
      multiMapping = multiMapping,
      mappingQuality = mappingQuality,
      threads = workers,
      outputFileName = countFileName
  }
  
  output {
    Array[File] trimStatsFiles = TrimGaloreSingle.statsFile
    File countFile = FeatureCounts.outputFile
    File countSummaryFile = FeatureCounts.summaryFile
  }
}