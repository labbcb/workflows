version 1.0
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/trimgalore/0.4.5/trimgalore.wdl" as trimgalore
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/bismark/0.19.0/bismark.wdl" as bismark

workflow WGBS {

    input {
        # Raw Sequencing Reads
        Array[File] filesR1
        Array[File] filesR2

        # Reference Genome Files
        Array[File] genomeFiles

        # Read Trimming
        Boolean illumina = true
        Int clipR1 = 8
        Int clipR2 = 8

        # Read Alignment
        Boolean unmapped = true
        Boolean ambiguous = true

        # Deduplicate Aligned Sequences
        Boolean paired = true
        Boolean bam = true

        # Methylation Quantification
        Boolean bedGraph = true
    }

    call bismark.GenomePreparation {
        input:
            genomeFiles = genomeFiles
    }

    scatter (pairedFiles in zip(filesR1, filesR2)) {
        call trimgalore.TrimGalorePaired {
            input:
                pairedFiles = pairedFiles,
                illumina = illumina,
                clipR1 = clipR1,
                clipR2 = clipR2
        }

        call bismark.BismarkPaired {
            input:
                genomeFiles = genomeFiles,
                indexFilesCT = GenomePreparation.indexFilesCT,
                indexFilesGA = GenomePreparation.indexFilesGA,
                pairedFiles = TrimGalorePaired.trimFiles,
                unmapped = unmapped,
                ambiguous = ambiguous
            }

        call bismark.Deduplicate {
            input:
                file = BismarkPaired.outputFile,
                paired = paired,
                bam = bam
        }

        call bismark.MethylationExtractor {
            input:
                file = Deduplicate.outputFile,
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
