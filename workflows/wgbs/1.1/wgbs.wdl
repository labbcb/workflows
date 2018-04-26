import "https://raw.githubusercontent.com/labbcb/rnnr/master/tools/trimgalore/0.4.4/trimgalore.wdl" as trimgalore
import "https://raw.githubusercontent.com/labbcb/rnnr/master/tools/bismark/0.19.0/bismark.wdl" as bismark

workflow WGBS {
  # Raw Sequencing Reads
  Array[File] filesR1
  Array[File] filesR2

  # Reference Genome Files
  Array[File] genomeFiles

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

  scatter (pairedFiles in zip(filesR1, filesR2)) {
 	  call trimgalore.TrimGalorePaired {
      input:
        pairedFiles = pairedFiles,
        trim1 = trim1,
        illumina = illumina
    }
  }

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
    Array[File] outputFile = MethylationExtractor.outputFile

    # Report files
    Array[Pair[File, File]] trimReportFiles = TrimGalorePaired.statsFiles
    Array[File] alignReportFile = BismarkPaired.reportFile
    Array[File] deduplicatedReportFiles = Deduplicate.reportFile

    # Intermediary Files
    #Array[File] indexFilesCT = GenomePreparation.indexFilesCT
    #Array[File] indexFilesGA = GenomePreparation.indexFilesGA
    #Array[File] alignFiles = BismarkPaired.outputFile
    #Array[Pair[File, File]] ambiguousFiles = BismarkPaired.ambiguousFiles
    #Array[Pair[File, File]] unmappedFiles = BismarkPaired.unmappedFiles
    #Array[File] deduplicatedFiles = Deduplicate.outputFile
  }
}
