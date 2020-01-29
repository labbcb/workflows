## Convert one or more BAM files to CRAM format and create index
##
## Requirements/expectations :
## - An array of input files in BAM format
## - Reference genome, in FASTA format, used to align the sequences
## Outputs :
## - An array of output files in CRAM format
## - An array of CRAM index files

workflow BamToCram {

    Array[File] array_bams
    File ref_fasta

    String? gitc_docker_override
    String gitc_docker = select_first([gitc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"])
    String? samtools_path_override
    String samtools_path = select_first([samtools_path_override, "samtools"])

    scatter (bam in array_bams) {
        call ConvertBamToCram {
            input:
                bam_file = bam,
                ref_fasta = ref_fasta,
                docker = gitc_docker,
                samtools_path = samtools_path
        }
    }

    output {
        Array[File] cram_files = ConvertBamToCram.cram_file
        Array[File] cram_index = ConvertBamToCram.cram_index
    }
}

task ConvertBamToCram {
    File bam_file
    File ref_fasta

    String docker
    String samtools_path
    Int? cpu

    String output_filename = basename(bam_file, ".bam") + ".cram"

    command <<<
        set -eo pipefail
        ${samtools_path} view -T ${ref_fasta} -C -o ${output_filename} ${bam_file}
        ${samtools_path} index ${output_filename}
    >>>

    runtime {
        docker: docker
        cpu: select_first([cpu, 1])
    }

    output {
        File cram_file = output_filename
        File cram_index = "${output_filename}.crai"
    }
}