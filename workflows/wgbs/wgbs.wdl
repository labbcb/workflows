import "https://raw.githubusercontent.com/labbcb/tool-rqc/master/rqc.wdl" as rqc
import "https://raw.githubusercontent.com/labbcb/tool-trimgalore/v0.4.4/trimgalore.wdl" as trimgalore
import "https://raw.githubusercontent.com/labbcb/tool-fastqc/v0.11.6/fastqc.wdl" as fastqc
import "https://raw.githubusercontent.com/labbcb/tool-bismark/v0.14.5/bismark.wdl" as bismark

workflow WGBS {

  # Raw Sequencing Reads
  Array[File] files
  Array[String] groups
  Array[File] filesR1
  Array[File] filesR2
  Array[String] groupsR1
  Array[String] groupsR2

  # Reference Genome Files
  Array[File] genomeFiles

  # Parallel Processing
  Int workers = 1
  
  # Quality Assessment
  String reportRawFilename = "reportRaw"
  String reportTrimFilename = "reportTrim"

  # Read Trimming
  Boolean trim1 = true
  Boolean illumina = true

  # Read Alignment
  Int seedmms = 1
  Boolean bowtie1 = true
  Boolean unmapped = true
  Boolean ambiguous = true

  # Deduplicate Aligned Sequences
  Boolean paired = true
  Boolean bam = true

  # Methylation Quantification
  Boolean bedGraph = true

  call rqc.Rqc as RqcRaw {
    input:
      files = files,
      groups = groups,
      workers = workers
   }
   
   scatter(file in files) {
        call fastqc.FastQC as FastQCRaw {
            input:
            	file = file
        }
    }

  scatter (pairedFiles in zip(filesR1, filesR2)) {
 	  call trimgalore.TrimGalorePaired {
      input:
        pairedFiles = pairedFiles,
        trim1 = trim1,
        illumina = illumina
    }
  }

#  call rqc.Rqc as RqcTrim {
#    input:
#      files = TrimGalorePaired.trimFiles
#      groups = groups,
#      workers = workers
#  }

  call bismark.GenomePreparation {
    input:
      genomeFiles = genomeFiles,
      bowtie1 = bowtie1
  }

  scatter (pairedFiles in TrimGalorePaired.trimFiles) {
    call bismark.BismarkPaired {
      input:
        genomeFiles = genomeFiles,
        indexFilesCT = GenomePreparation.indexFilesCT,
        indexFilesGA = GenomePreparation.indexFilesGA,
        pairedFiles = pairedFiles,
        seedmms = seedmms,
        bowtie1 = bowtie1,
        unmapped = unmapped,
        ambiguous = ambiguous
      }
  }

  scatter (file in BismarkPaired.outputFile) {
    call bismark.Deduplicate {
      input:
        file = file,
        paired = paired,
        bam = bam
    }
  }

  scatter (file in Deduplicate.outputFile) {
    call bismark.MethylationExtractor {
      input:
        file = file,
        paired = paired,
        bedGraph = bedGraph
    }
  }

  output {
    File integratedReportRaw = RqcRaw.report
    Array[File] reportsRaw = FastQCRaw.report
    Array[Pair[File, File]] trimStatsFiles = TrimGalorePaired.statsFiles
    Array[File] alignFiles = BismarkPaired.outputFile
    Array[File] reportFiles = BismarkPaired.reportFile
    Array[File] deduplicatedReportFiles = Deduplicate.reportFile
    Array[File] outputFile = MethylationExtractor.outputFile
    
    # Intermediary Files
    #Array[File] indexFilesCT = GenomePreparation.indexFilesCT
    #Array[File] indexFilesGA = GenomePreparation.indexFilesGA
    #Array[Pair[File, File]] ambiguousFiles = BismarkPaired.ambiguousFiles
    #Array[Pair[File, File]] unmappedFiles = BismarkPaired.unmappedFiles
    #Array[File] deduplicatedFiles = Deduplicate.outputFile
  }
}
