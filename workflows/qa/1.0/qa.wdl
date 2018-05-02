import "tools/rqc/1.14.0/rqc.wdl" as rqc
import "tools/fastqc/0.11.7/fastqc.wdl" as fastqc

workflow QA {
    Array[File]+ files
    Array[String]+ groups

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
        File rqcReportFile = Rqc.report
        Array[File] fastqcReportFiles = FastQC.report
    }
}
