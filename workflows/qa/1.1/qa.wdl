version 1.0

import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/rqc/1.16.2/rqc.wdl" as rqc
import "https://raw.githubusercontent.com/labbcb/workflows/master/tools/fastqc/0.11.8/fastqc.wdl" as fastqc

workflow QA {

    input {
        Array[File]+ files
        Array[String]+ groups
    }

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
