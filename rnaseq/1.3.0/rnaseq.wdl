version 1.0

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
        call TrimGalorePaired {
            input:
                pairedFiles = (filesR1[idx], filesR2[idx]),
                gzip = gzip
        }

        call AlignReads {
            input:
                indexFiles = indexFiles,
                pairedFiles = TrimGalorePaired.trimFiles
        }

        call Count {
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

task TrimGaloreSingle {

    input {
        File file
        String? adapter
        Int? maxLength
        Int? maxN
        Int? clipR1
        Int? threePrimeClipR1

        # General options
        Int quality = 20
        Boolean phred33 = true
        Boolean phred64 = false
        Boolean illumina = false
        Boolean nextera = false
        Boolean smallRNA = false
        Int stringency = 1
        Float errorRate = 0.1
        Boolean gzip = true
        Int length = 20
        Boolean trimN = false

        # RRBS-specific options (MspI digested material)
        Boolean rrbs = false
        Boolean nonDirectional = false
        Boolean keep = false
	}

	command {
		trim_galore ${file} \
			--quality ${quality} \
			${true='--phred33' false='' phred33} \
			${true='--phred64' false='' phred64} \
			${'--adapter ' + adapter} \
			${true='--illumina' false='' illumina} \
			${true='--nextera' false='' nextera} \
			${true='--small_rna' false='' smallRNA} \
			${'--max_length ' + maxLength} \
			--stringency ${stringency} \
			-e ${errorRate} \
			${true='--gzip' false='--dont_gzip' gzip} \
			--length ${length} \
			${'--max_n ' + maxN} \
			${true='--trim_n' false='' trimN} \
			${'--clip_R1 ' + clipR1} \
			${'--three_prime_clip_R1 ' + threePrimeClipR1} \
			${true='--rrbs' false='' rrbs} \
			${true='--non_directional' false='' nonDirectional} \
			${true='--keep' false='' keep} \
	}

	output {
		File trimFile = sub(basename(file), ".fastq(.gz)?", "_trimmed.fq") + if (gzip) then ".gz" else ""
		File statsFile = basename(file) + "_trimming_report.txt"
	}

	runtime {
		docker: "welliton/trimgalore:0.5.0"
	}
}

task AlignReads {

    input {
        Array[File] indexFiles
        Pair[File, File] pairedFiles
    }

    command {
        mkdir index
        ln -s ${sep=' ' indexFiles} -t index
        STAR --runMode alignReads --genomeDir index \
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