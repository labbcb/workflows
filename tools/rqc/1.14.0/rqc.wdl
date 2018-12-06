version 1.0

task Rqc {

    input {
        Array[File] files
        Array[String] groups
        Boolean sample = true
        Int reads = 1000000
        Int workers = 1
        String reportFile = "rqc_report"
        String rdsFile = "rqc.rds"
    }

	command {
		rqc --files ${sep=',' files} \
			--groups ${sep=',' groups} \
			${true='--sample' false='' sample} \
			--reads ${reads} \
			--workers ${workers} \
			--report_file ${reportFile} \
			--rds_file ${rdsFile}
	}

	output {
		File report = "${reportFile}.html"
		File rds = "${rdsFile}"
	}

	runtime {
		docker: "welliton/rqc:1.14.0"
	}

	parameter_meta {
		files: "List of sequencing files in FASTQ or BAM format."
		groups: "List of sample groups."
		sample: "Process a random sample from files."
		reads: "Number of sequences to read from each input file. This represents sample size --sample flag is present. Default value is 1000000 reads."
		workers: "Number of processing units for parallelization."
		reportFile: "Report file name without extension."
		rdsFile: "RDS file name with extension."
	}

	meta {
		author: "Welliton Souza"
		email: "well309@gmail.com"
	}
}
