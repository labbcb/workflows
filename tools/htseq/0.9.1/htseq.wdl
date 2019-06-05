version 1.0

task Count {

    input {
        File file
        File gtfFile
        String format = "sam"
        String stranded = "yes"
        String destination = "stdout"
    }

    command {
        htseq-count --format ~{format} --stranded ~{stranded} \
            ~{file} ~{gtfFile} > ~{destination}
    }

    output {
        File countFile = destination
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
        String destination = "stdout"
    }

    command {
        htseq-count --format ~{format} --stranded ~{stranded} \
            ~{sep=' ' files} > ~{destination}
    }

    output {
        File countFile = destination
    }

    runtime {
        docker: "welliton/htseq:0.9.1"
    }
}
