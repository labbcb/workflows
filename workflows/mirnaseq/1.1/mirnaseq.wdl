version 1.0
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/trimgalore/0.5.0/trimgalore.wdl" as trimgalore
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/bowtie/1.2.2/bowtie.wdl" as bowtie
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/subread/1.6.2/subread.wdl" as subread

workflow MirnaSeq {

    input {
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
      Int seedmms = 2
      Boolean all = true
      Boolean best = true
      Boolean strata = true
      Boolean sam = true
      String featureType = "miRNA"
      String attributeType = "Alias"
      Boolean overlaps = true
      Boolean multiMapping = true
      Int mappingQuality = 20
    }
  
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
        maxN = maxN
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
    File countFile = FeatureCounts.outputFile
    Array[File] trimStatsFiles = TrimGaloreSingle.statsFile
    Array[File] alignStatsFiles = Align.statsFile
    File countSummaryFile = FeatureCounts.summaryFile
  }
}