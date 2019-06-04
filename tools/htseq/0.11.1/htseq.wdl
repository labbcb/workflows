version 1.0

task Count {

    input {
        File file
        File gtfFile
        String destination
        String format = "sam"
        String stranded = "yes"
    }

    command {
        htseq-count \
            --format ~{format} \
            --stranded ~{stranded} \
            ~{file} ~{gtfFile} > ~{destination}
    }

    output {
        File countFile = destination
    }

    runtime {
        docker: "welliton/htseq:0.11.1"
    }
}

task CountMultipleFiles {

    input {
        Array[File]+ files
        File gtfFile
        String destination
        String format = "sam"
        String stranded = "yes"
    }

    command {
        htseq-count \
            --format ~{format} \
            --stranded ~{stranded} \
            ~{sep=' ' files} ~{gtfFile} > ~{destination}
    }

    output {
        File countFile = destination
    }

    runtime {
        docker: "welliton/htseq:0.11.1"
    }
}
