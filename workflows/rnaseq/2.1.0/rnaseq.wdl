# Copyright BCBLab, 2019
# RNA-seq data processing for differential gene expression analysis
# Inputs
# - Arrays of paths to aired-end FASTQ files (fastq_1, fastq_2)
# - Array of sample names (samples_names)
# - Paths of genome file (genome), and STAR index files (sa, sa_index, chr_length, chr_name, chr_name_length, chr_start, genome_parameters)
# - Path to GTF file (gtf), must be the same file used in generating genome index
# Outputs
# - Array of raw count files
version 1.0

workflow RNAseq {

    input {
        Array[File] fastq_1
        Array[File] fastq_2

        Array[File] sample_names

        File genome
        File sa
        File sa_index
        File chr_length
        File chr_name
        File chr_name_length
        File chr_start
        File genome_parameters

        File gtf
        String stranded = "reverse"
    }

    scatter (idx in range(length(sample_names))) {

        String sample_name = sample_names[idx]

        call Trim {
            input:
                fastq_1 = fastq_1[idx],
                fastq_2 = fastq_2[idx]
        }

        call Align {
            input:
                fastq_1 = Trim.trim_1,
                fastq_2 = Trim.trim_2,
                sample_name = sample_name,
                Genome = genome,
                SA = sa,
                SAindex = sa_index,
                chrLength = chr_length,
                chrName = chr_name,
                chrNameLength = chr_name_length,
                chrStart = chr_start,
                genomeParameters = genome_parameters
        }

        call Count {
            input:
                bam_file = Align.align,
                gtf_file = gtf,
                stranded = stranded,
                destination = sample_name + ".txt"
        }
    }

    output {
        Array[File] counts = Count.count
    }

}

task Trim {

    input {
        File fastq_1
        File fastq_2
        Boolean illumina = true
        Boolean gzip = true
	}

	command {
		trim_galore ~{fastq_1} ~{fastq_2} --illumina --gzip --paired
	}

	output {
		File trim_1 = basename(fastq_1, ".fastq.gz") + "_val_1.fq.gz"
		File trim_2 = basename(fastq_2, ".fastq.gz") + "_val_2.fq.gz"
		File stats_1 = basename(fastq_1) + "_trimming_report.txt"
        	File stats_2 = basename(fastq_2) + "_trimming_report.txt"
	}

	runtime {
		docker: "welliton/trimgalore:0.6.4"
	}
}

task Align {

    input {
        File fastq_1
        File fastq_2
        String sample_name

        File Genome
        File SA
        File SAindex
        File chrLength
        File chrName
        File chrNameLength
        File chrStart
        File genomeParameters
    }

    command {
        STAR --runMode alignReads \
            --genomeDir ~{sub(Genome, basename(Genome), "")} \
            --readFilesIn ~{fastq_1} ~{fastq_2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ~{sample_name + "_"}
    }

    output {
        File align = "${sample_name}_Aligned.sortedByCoord.out.bam"
    }

    runtime {
        docker: "welliton/star:2.7.3a"
    }
}

task Count {

    input {
        File bam_file
        File gtf_file
        String stranded
        String destination
    }

    command {
        htseq-count --format bam --stranded ~{stranded} \
            ~{bam_file} ~{gtf_file} > ~{destination}
    }

    output {
        File count = destination
    }

    runtime {
        docker: "welliton/htseq:0.11.1"
    }
}
