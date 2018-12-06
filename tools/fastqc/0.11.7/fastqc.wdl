version 1.0

task FastQC {

    input {
        File file
        Boolean? casava
        Boolean? nano
        Boolean? nofilter
        Boolean? extract = false
        Boolean? nogroup
        Int? min_length
        String? format
        Int? threads
        File? contaminants
        File? adapters
        File? limits
        Int? kmers = 7
        Boolean? quiet
    }

    command {
        fastqc --outdir "." \
            ${true='--casava' false='' casava} \
            ${true='--nano' false='' nano} \
            ${true='--nofilter' false='' nofilter} \
            ${true='--extract' false='--noextract' extract} \
            ${true='--nogroup' false='' nogroup} \
            ${'--min_length ' + min_length} \
            ${'--format ' + format} \
            ${'--contaminants ' + contaminants} \
            ${'--adapters ' + adapters} \
            ${'--limits ' + limits} \
            ${'--kmers ' + kmers} \
            ${true='--quiet' false='' quiet} \
            ${file}
    }

    runtime {
        docker: "welliton/fastqc:0.11.7"
    }

    output {
        File report = sub(basename(file), "\\..+", "_fastqc.html")
        File reportZip = sub(basename(file), "\\..+", "_fastqc.zip")
    }

    parameter_meta {
        file: "Input file."
        casava: "Files come from raw casava output."
        nano: "Files come from naopore sequences and are in fast5 format."
        nofilter: "If running with --casava then don't remove read flagged by casava as poor quality when performing the QC analysis."
        extract: "If set then the zipped output file will be uncompressed in the same directory after it has been created."
        nogroup: "Disable grouping of bases for reads >50bp."
        min_length: "Sets an artificial lower limit on the length of the sequence to be shown in the report."
        format: "Bypasses the normal sequence file format detection and forces the program to use the specified format. Valid formats are bam, sam, bam_mapped, sam_mapped and fastq."
        contaminants: "Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against."
        adapters: "Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library."
        limits: "Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules."
        kmers: "Specifies the length of Kmer to look for in the Kmer content module. Specified Kmer length must be between 2 and 10. Default length is 7 if not specified."
        quiet: "Supress all progress messages on stdout and only report errors."
    }

    meta {
        author: "Welliton Souza"
        email: "well309@gmail.com"
    }
}
