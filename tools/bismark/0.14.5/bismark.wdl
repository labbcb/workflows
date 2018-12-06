version 1.0

task GenomePreparation {

    input {
      Array[File] genomeFiles
      Boolean bowtie1 = false
      Boolean singleFasta = false
      Boolean genomicComposition = false
    }

  command {
    mv ${sep=' ' genomeFiles} .
    bismark_genome_preparation . \
      ${true='--bowtie1' false='--bowtie2' bowtie1} \
      ${true='--single_fasta' false='' singleFasta} \
      ${true='--genomic_composition' false='' genomicComposition}
  }

  output {

    Array[File] indexFilesCT = glob("Bisulfite_Genome/CT_conversion/*")
    Array[File] indexFilesGA = glob("Bisulfite_Genome/GA_conversion/*")
  }

  runtime {
		docker: "welliton/bismark:0.14.5"
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
      Int seedmms = 1
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
      ${true='--bowtie1' false='' bowtie1} \
      --seedmms ${seedmms} \
      ${true='--unmapped' false='' unmapped} \
      ${true='--ambiguous' false='' ambiguous}
  }

  output {
    File outputFile = "${basename(pairedFiles.left)}_bismark_pe.bam"
    File reportFile = "${basename(pairedFiles.left)}_bismark_PE_report.txt"
    Pair[File, File] ambiguousFiles = (
      "${basename(pairedFiles.left)}_ambiguous_reads_1.fq.gz",
      "${basename(pairedFiles.right)}_ambiguous_reads_2.fq.gz")
    Pair[File, File] unmappedFiles = (
      "${basename(pairedFiles.left)}_unmapped_reads_1.fq.gz",
      "${basename(pairedFiles.right)}_unmapped_reads_2.fq.gz")
  }

  runtime {
    docker: "welliton/bismark:v0.14.5"
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
    docker: "welliton/bismark:0.14.5"
  }
}

task MethylationExtractor {

    input {
      File file
      Boolean paired = false
      Boolean bedGraph = false
    }

  command {
    bismark_methylation_extractor ${file} \
      ${true='--paired-end' false='--single-end' paired} \
      ${true='--bedGraph' false='' bedGraph}
  }

  output {
    File outputFile = sub(basename(file), ".bam$", ".bismark.cov.gz")
    File bedGraphFile = sub(basename(file), ".bam$", ".bedGraph.gz")
  }

  runtime {
    docker: "welliton/bismark:0.14.5"
  }
}
