version 1.0
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/trimgalore/0.5.0/trimgalore.wdl" as trimgalore
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/star/2.5.3a/star.wdl" as star
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/htseq/0.11.1/htseq.wdl" as htseq

workflow RNAseq {

    input {
        # Sample names
        Array[String] sampleNames

        # Raw Sequencing Reads
        Array[File] filesR1
        Array[File] filesR2

        # Reference Genome Files
        Array[File] indexFiles
        File gtfFile

        # Read Trimming
        Boolean illumina = true
        Boolean gzip = false

        # Gene Quantification
        String format = "bam"
        String stranded = "reverse"
    }


    scatter (idx in range(length(sampleNames))) {
        call trimgalore.TrimGalorePaired {
            input:
                pairedFiles = (filesR1[idx], filesR2[idx]),
                gzip = gzip
        }

        call star.AlignReads {
            input:
                indexFiles = indexFiles,
                pairedFiles = TrimGalorePaired.trimFiles
        }

        call htseq.Count {
            input:
                file = AlignReads.alignFile,
                destination = sampleNames[idx] + ".txt",
                gtfFile = gtfFile,
                format = format,
                stranded = stranded
        }
    }

    output {
        Array[File] countFiles = Count.countFile
    }
}