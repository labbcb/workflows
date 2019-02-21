version 1.0

task Count {

    input {
        File file
        File gtfFile
        String format = "sam"
        String stranded = "yes"
    }

    command {
        htseq-count --format ${format} --stranded ${stranded} \
            ${file} ${gtfFile}
    }

    output {
        File countFile = stdout()
    }

    runtime {
        docker: "welliton/htseq:0.9.1"
    }
}

task CountMultipleFiles {

    input {
        Array[File]+ files
        File gtfFile
        String format = "sam"
        String stranded = "yes"
    }

    command {
        htseq-count --format ${format} --stranded ${stranded} \
            ${sep=' ' files} ${gtfFile}
    }

    output {
        File countFile = stdout()
    }

    runtime {
        docker: "welliton/htseq:0.9.1"
    }
}
