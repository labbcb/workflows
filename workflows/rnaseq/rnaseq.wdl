import "http://raw.githubusercontent.com/labbcb/tool-rqc/master/rqc.wdl" as rqc
import "http://raw.githubusercontent.com/labbcb/tool-fastqc/master/v0.11.6/fastqc.wdl" as fastqc

workflow RNAseq {
	Array[File] files

	call rqc.Rqc {
		input:
			files = files
	}

	scatter(file in files) {
    	call fastqc.FastQC as fastqc_raw {
        	input:
				file = file
    	}
	}

	output {
		File qa_report_raw = Rqc.qc_report
		Array[File] fastqc_raw_reports = fastqc_raw.report
	}
}