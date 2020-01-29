version 1.0

workflow MirnaSeq {

    input {
        Array[File] files
        Array[File] indexFiles
        String ebwtBase
        File annotationFile

        String countFileName = "counts.txt"
        Int workers = 1
        Boolean smallRNA = true
        Int length = 18
        Int maxN = 2
        Int maxLength = 30
        Int seedmms = 2
        Boolean all = true
        Boolean best = true
        Boolean strata = true
        Boolean sam = true
        String featureType = "miRNA"
        String attributeType = "Alias"
        Boolean overlaps = true
        Boolean multiMapping = true
        Int mappingQuality = 20
    }
    
    scatter (file in files) {
        call TrimGaloreSingle {
            input:
                file = file,
                smallRNA = smallRNA,
                length = length,
                maxLength = maxLength,
                maxN = maxN
        }
        call Align {
            input:
                indexFiles = indexFiles,
                ebwtBase = ebwtBase,
                file = TrimGaloreSingle.trimFile,
                outputFileName = sub(basename(file), "_R1.fastq.gz", ".sam"),
                seedmms = seedmms,
                all = all,
                best = best,
                strata = strata,
                sam = sam,
                threads = workers
        }
    }
    
    call FeatureCounts {
        input:
            files = Align.outputFile,
            annotationFile = annotationFile,
            featureType = featureType,
            attributeType = attributeType,
            overlaps = overlaps,
            multiMapping = multiMapping,
            mappingQuality = mappingQuality,
            threads = workers,
            outputFileName = countFileName
    }
    
    output {
        File countFile = FeatureCounts.outputFile
        Array[File] trimStatsFiles = TrimGaloreSingle.statsFile
        Array[File] alignStatsFiles = Align.statsFile
        File countSummaryFile = FeatureCounts.summaryFile
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


task Align {

    input {
        Array[File] indexFiles
        String ebwtBase
        File file
        String outputFileName
        Int seedmms = 2
        Boolean all = false
        Boolean best = false
        Boolean strata = false
        Boolean sam = false
        Int threads = 1
    }

    command {
        ln -s ${sep=' ' indexFiles} -t .
        bowtie \
            --seedmms ${seedmms} \
            ${true='--all' false='' all} \
            ${true='--best' false='' best} \
            ${true='--strata' false='' strata} \
            ${true='--sam' false='' strata} \
            --threads ${threads} \
            ${ebwtBase} \
            ${file} \
            ${outputFileName}
    }
    
    output {
        File outputFile = outputFileName
        File statsFile = stderr()
    }
    
    runtime {
        docker: "welliton/bowtie:1.2.2"
    }
    
    parameter_meta {
        indexFiles: "Index files (*.ebwt) created by BuidIndex task."
        ebwtBase: "The basename of the index to be searched."
        file: "A single file containing unpaired reads to be aligned."
        outputFileName: "File to write alignments to."
        seedmms: "Maximum number of mismatches permitted in the seed. This may be 0, 1, 2 or 3 and the default is 2."
        all: "Report all valid alignments per read or pair."
        best: "Make Bowtie guarantee that reported singleton alignments are best in terms of stratum."
        strata: "If many valid alignments exist and are reportable and they fall into more than one alignment stratum, report only those alignments that fall into the best stratum."
        sam: "Print alignments in SAM format."
        threads: "Launch parallel search threads (default: 1)."
    }
    
    meta {
        author: "Welliton Souza"
        email: "well309@gmail.com"
    }
}

task FeatureCounts {

    input {
        File annotationFile
        Array[File] files
        String outputFileName

        String format = "GTF"
        String featureType = "exon"
        String attributeType = "gene_id"
        Boolean featureLevel = false
        Boolean overlaps = false
        Boolean multiMapping = false
        Int mappingQuality = 0
        Int threads = 1
    }

    command {
        featureCounts \
            -F ${format} \
            -t ${featureType} \
            -g ${attributeType} \
            ${true='-f' false='' featureLevel} \
            ${true='-O' false='' overlaps} \
            ${true='-M' false='' multiMapping} \
            -Q ${mappingQuality} \
            -T ${threads} \
            -a ${annotationFile} \
            -o ${outputFileName} \
            ${sep=' ' files}
    }

    output {
        File outputFile = outputFileName
        File summaryFile = "${outputFileName}.summary"
    }

    runtime {
        docker: "welliton/subread:1.6.2"
    }
}