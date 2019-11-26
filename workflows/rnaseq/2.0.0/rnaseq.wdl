## Copyright BCBLab, 2019
## RNA-seq data processing for differential gene expression analysis
version 1.0

workflow RNAseq {

    input {
        Array[File] fastq_1
        Array[File] fastq_2

        String suffix_1 = "_R1.fastq.gz"

        File Genome
        File SA
        File SAindex
        File chrLength
        File chrName
        File chrNameLength
        File chrStart
        File genomeParameters

        File gtf
        String stranded = "reverse"
    }

    scatter (idx in range(length(files_1))) {

        String sample_name = basename(files_1[idx], suffix_1) + ".txt"

        call Trim {
            input:
                fasfq_1 = fastq_1[idx],
                fastq_2 = fastq_2[idx]
        }

        call Align {
            input:
                fastq_1 = Trim.trim_1,
                fastq_2 = Trim.trim_2,
                sample_name = sample_name,
                Genome = Genome,
                SA = SA,
                SAindex = SAindex,
                chrLength = chrLength,
                chrName = chrName,
                chrNameLength = chrNameLength,
                chrStart = chrStart,
                genomeParameters = genomeParameters,
                gtf = gtf
        }

        call Count {
            input:
                file = Align.align,
                gtf = gtf,
                format = format,
                stranded = stranded,
                destination = sample_name
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
		File trim_2	= basename(fastq_2, ".fastq.gz") + "_val_2.fq"
		File stats_1 = basename(fastq_1) + "_trimming_report.txt"
        File stats_2 = basename(fastq_2) + "_trimming_report.txt"
	}

	runtime {
		docker: "welliton/trimgalore:0.5.0"
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

        File gtf
    }

    command {
        STAR --runMode alignReads \
            --genomeDir ~{sub(Genome, basename(Genome), "")} \
            --readFilesIn ~{fastq_1} ~{fastq_2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbGTFfile ~{gtf} \
            --outFileNamePrefix ~{sample_name + "_"}
    }

    output {
        File align = "${sample_name}_Aligned.sortedByCoord.out.bam"
    }

    runtime {
        docker: "welliton/star:2.5.3a"
    }
}

task Count {

    input {
        File file
        File gtf
        String format
        String stranded
        String destination
    }

    command {
        htseq-count --format bam --stranded ~{stranded} \
            ~{file} ~{gtf} > ~{destination}
    }

    output {
        File count = destination
    }

    runtime {
        docker: "welliton/htseq:0.9.1"
    }
}