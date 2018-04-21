task TrimGaloreSingle {
	# General options
	File file
	Int quality = 20
	Boolean phred33 = true
	Boolean phred64 = false
	String? adapter
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
	Int? threePrimeClipR1
	
	# RRBS-specific options (MspI digested material)
	Boolean rrbs = false
	Boolean nonDirectional = false
	Boolean keep = false
	
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
		docker: "welliton/trimgalore:v0.4.4"
	}
}

task TrimGalorePaired {
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
		docker: "welliton/trimgalore:v0.4.4"
	}
}