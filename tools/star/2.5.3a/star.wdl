version 1.0

task GenomeGenerate {

    input {
        Array[File] genomeFiles
        File? gtfFile
    }

    command {
        STAR --runMode genomeGenerate --genomeDir . \
            --genomeFastaFiles ${sep=' ' genomeFiles} \
            ${'--sjdbGTFfile ' + gtfFile}
    }

    output {
        Array[File] indexFiles = ["Genome", "Log.out", "SA", "SAindex", "chrLength.txt", "chrName.txt", "chrNameLength.txt", "chrStart.txt", "genomeParameters.txt"]
    }

    runtime {
        docker: "welliton/star:2.5.3a"
    }
}

task AlignReads {

    input {
        Array[File] indexFiles
        Pair[File, File] pairedFiles
    }

    command {
        ln -s ${sep=' ' indexFiles} -t .
        STAR --runMode alignReads --genomeDir . \
            --readFilesIn ${pairedFiles.left} ${pairedFiles.right} \
            --outSAMtype BAM SortedByCoordinate
    }

    output {
        File alignFile = "Aligned.sortedByCoord.out.bam"
        File logFinalFile = "Log.final.out"
        File logFile = "Log.out"
        File logProgressFile = "Log.progress.out"
        File sjFile = "SJ.out.tab"
    }

    runtime {
        docker: "welliton/star:2.5.3a"
    }
}
