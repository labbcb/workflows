version 1.0

task Reaper {

    input {
        File inputFile
        String geometry # no-bc (no barcode), 5p-bc (5' barcode) or 3p-bc (3' barcode)
        String? adapter3p
        String? tabu
        Int minLength = 0
        Boolean nozip = false
        String outputBaseName = basename(inputFile)
    }

    command {
        reaper \
            -i ${inputFile} \
            -geom ${geometry} \
            ${'-3pa ' + adapter3p} \
            ${'-tabu ' + tabu} \
            -clean-length ${minLength} \
            ${true='--nozip' false='' nozip} \
            ${'-basename ' + outputBaseName}
    }

    output {
        File outputFile = "${outputBaseName}.lane.clean.gz"
    }

    runtime {
        docker: "welliton/kraken:13-274"
    }
}
