import "https://raw.githubusercontent.com/labbcb/tool-rqc/master/rqc.wdl" as rqc
import "https://raw.githubusercontent.com/labbcb/tool-fastqc/master/v0.11.6/fastqc.wdl" as fastqc

workflow QA {
    Array[File]+ files
    Array[String] groups

    call rqc.Rqc {
    	input:
    		files = files,
    		groups = groups
    }

    scatter(file in files) {
        call fastqc.FastQC {
            input:
            	file = file
        }
    }

    output {
        File integratedReport = Rqc.report
        Array[File] reports = FastQC.report
    }
}
