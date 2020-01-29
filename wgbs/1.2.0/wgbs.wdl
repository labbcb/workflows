version 1.0

workflow WGBS {

    input {
        # Raw Sequencing Reads
        Array[File] filesR1
        Array[File] filesR2

        # Reference Genome Files
        Array[File] genomeFiles
        Array[File] indexFilesCT
        Array[File] indexFilesGA

        # Read Trimming
        Boolean illumina = true
        Int clipR1 = 8
        Int clipR2 = 8

        # Read Alignment
        Boolean unmapped = true
        Boolean ambiguous = true

        # Deduplicate Aligned Sequences
        Boolean paired = true
        Boolean bam = true

        # Methylation Quantification
        Boolean bedGraph = true
    }

    scatter (pairedFiles in zip(filesR1, filesR2)) {
        call TrimGalorePaired {
            input:
                pairedFiles = pairedFiles,
                illumina = illumina,
                clipR1 = clipR1,
                clipR2 = clipR2
        }

        call BismarkPaired {
            input:
                genomeFiles = genomeFiles,
                indexFilesCT = indexFilesCT,
                indexFilesGA = indexFilesGA,
                pairedFiles = TrimGalorePaired.trimFiles,
                unmapped = unmapped,
                ambiguous = ambiguous
            }

        call Deduplicate {
            input:
                file = BismarkPaired.outputFile,
                paired = paired,
                bam = bam
        }

        call MethylationExtractor {
            input:
                file = Deduplicate.outputFile,
                paired = paired,
                bedGraph = bedGraph
        }
    }

    output {
        Array[File] outputFile = MethylationExtractor.outputFile

        # Report files
        Array[Pair[File, File]] trimReportFiles = TrimGalorePaired.statsFiles
        Array[File] alignReportFile = BismarkPaired.reportFile
        Array[File] deduplicatedReportFiles = Deduplicate.reportFile

        # Intermediary Files
        #Array[File] alignFiles = BismarkPaired.outputFile
        #Array[Pair[File, File]] ambiguousFiles = BismarkPaired.ambiguousFiles
        #Array[Pair[File, File]] unmappedFiles = BismarkPaired.unmappedFiles
        #Array[File] deduplicatedFiles = Deduplicate.outputFile
    }
}

task TrimGalorePaired {

    input {
        # General options
        Pair[File, File] pairedFiles
        Int quality = 20
        Boolean phred33 = true
        Boolean phred64 = false
        String? adapter
        String? adapter2
        Boolean illumina = false
        Boolean nextera = false
        Boolean smallRNA = false
        Int? maxLength
        Int stringency = 1
        Float errorRate = 0.1
        Boolean gzip = true
        Int length = 20
        Int? maxN
        Boolean trimN = false
        Int? clipR1
        Int? clipR2
        Int? threePrimeClipR1
        Int? threePrimeClipR2
	}

	# RRBS-specific options (MspI digested material)
	Boolean rrbs = false
	Boolean nonDirectional = false
	Boolean keep = false

	# Paired-end specific options
	Boolean trim1 = false

	command {
		trim_galore ${pairedFiles.left} ${pairedFiles.right} \
			--quality ${quality} \
			${true='--phred33' false='' phred33} \
			${true='--phred64' false='' phred64} \
			${'--adapter ' + adapter} \
			${'--adapter2 ' + adapter2} \
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
			${'--clip_R2 ' + clipR2} \
			${'--three_prime_clip_R1 ' + threePrimeClipR1} \
			${'--three_prime_clip_R2 ' + threePrimeClipR2} \
			${true='--rrbs' false='' rrbs} \
			${true='--non_directional' false='' nonDirectional} \
			${true='--keep' false='' keep} \
			--paired \
			${true='--trim1' false='' trim1}
	}

	output {
		Pair[File, File] trimFiles = (
			sub(basename(pairedFiles.left), ".fastq(.gz)?", "_val_1.fq") + if (gzip) then ".gz" else "",
			sub(basename(pairedFiles.right), ".fastq(.gz)?", "_val_2.fq") + if (gzip) then ".gz" else "")
		Pair[File, File] statsFiles = (
			basename(pairedFiles.left) + "_trimming_report.txt",
			basename(pairedFiles.right) + "_trimming_report.txt")
	}

	runtime {
		docker: "welliton/trimgalore:0.4.5"
	}
}

task BismarkPaired {

    input {
        Array[File] genomeFiles
        Array[File] indexFilesCT
        Array[File] indexFilesGA
        Pair[File, File] pairedFiles

        Boolean fastq = true

        # Alignment
        Boolean bowtie1 = false
        Int? seedmms

        # Output
        Boolean nonDirectional = false
        Boolean unmapped = false
        Boolean ambiguous = false
    }

    command {
        mv ${sep=' ' genomeFiles} .
        mkdir Bisulfite_Genome
        mkdir Bisulfite_Genome/CT_conversion/
        mkdir Bisulfite_Genome/GA_conversion/
        mv ${sep=' ' indexFilesCT} Bisulfite_Genome/CT_conversion/
        mv ${sep=' ' indexFilesGA} Bisulfite_Genome/GA_conversion/

        bismark . \
            -1 ${pairedFiles.left} -2 ${pairedFiles.right} \
            ${true='--fastq' false='--fasta' fastq} \
            ${true='--bowtie1' false='--bowtie2' bowtie1} \
            ${'--seedmms ' + seedmms} \
            ${true='--non_directional' false='' nonDirectional} \
            ${true='--unmapped' false='' unmapped} \
            ${true='--ambiguous' false='' ambiguous}
    }

    output {
        File outputFile = glob("*_pe.bam")[0]
        File reportFile = glob("*_PE_report.txt")[0]
        Pair[File, File] ambiguousFiles = (
            "${basename(pairedFiles.left)}_ambiguous_reads_1.fq.gz",
            "${basename(pairedFiles.right)}_ambiguous_reads_2.fq.gz")
        Pair[File, File] unmappedFiles = (
            "${basename(pairedFiles.left)}_unmapped_reads_1.fq.gz",
            "${basename(pairedFiles.right)}_unmapped_reads_2.fq.gz")
    }

    runtime {
        docker: "welliton/bismark:0.19.0"
    }
}

task Deduplicate {

    input {
        File file
        Boolean paired = true
        Boolean bam = false
    }

    command {
        mv ${file} .
        deduplicate_bismark ${basename(file)} \
            ${true='--paired' false='--single' paired} \
            ${true='--bam' false='' bam}
    }

    output {
        File outputFile = sub(basename(file), ".bam$", ".deduplicated.bam")
        File reportFile = sub(basename(file), ".bam$", ".deduplication_report.txt")
    }

    runtime {
        docker: "welliton/bismark:0.19.0"
    }
}

task MethylationExtractor {

    input {
        File file
        Boolean paired = false
        Boolean bedGraph = false
        Boolean counts = false
        Boolean zeroBased = false
    }

    command {
        bismark_methylation_extractor ${file} \
            ${true='--paired-end' false='--single-end' paired} \
            ${true='--bedGraph' false='' bedGraph} \
            ${true='--counts' false='' counts} \
            ${true='--zero_based' false='' zeroBased}
    }

    output {
        File outputFile = sub(basename(file), ".bam$", ".bismark.cov.gz")
        File bedGraphFile = sub(basename(file), ".bam$", ".bedGraph.gz")
    }

    runtime {
        docker: "welliton/bismark:0.19.0"
    }
}