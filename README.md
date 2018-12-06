# Reproducible Workflows for Genomics Research

This repository contains file descriptions for common bioinformatics tools and workflows.
These workflows have been developed and used by the Biostatistics and Computational Laboratory (BCBLab).
Descriptions of tools and workflows were written using reproducible languages for reproducible research:
Common Workflow Language (CWL) and Workflow Description Language (WDL).
Newest versions are mainly written in WDL.

## Tools

Dockerfiles and CWL/WDL files for tools are in `tools` directory.
Each subdirectory for one tool.
Inside each subdirectory there are many directories representing tool version.
The version must be equal to tool version downloaded by Dockerfile.
Dockerfiles are tagged according to tool version.
Each version must have its own Dockerfile and CWL/WDL files.
CWL/WDL must use the right version of Docker image according to tool version.

| Tool          | Versions |
| ------------- | -------- |
| bcl2fastq     | 1.8.4, 2.20.0 |
| Bioconductor  | 3.6, 3.7 |
| Bismark       | 0.14.5, 0.19.0 |
| Bowtie        | 1.2.0, 1.2.2 | 
| BWA           | 0.7.15, |
| Cutadapt      | 1.13, 1.16 |
| FastQC        | 0.11.6 0.11.7 |
| GATK          | 3.6-0 |
| HTSeq count   | 0.9.1 |
| Reaper        | 13-274 |
| Picard        | 2.6.0 |
| R             | 3.4.2, 3.5.0 |
| Rqc           | 1.14.0 |
| Samtools      | 1.5 |
| Star          | 2.5.3a |
| featureCounts | 1.5.2, 1.6.2 |
| TrimGalore!   | 0.4.4, 0.4.5, 0.5.0 |


## Workflows

Workflows import tool descriptions allowing reuse of these files.

| Workflow  | Versions | Description |
| --------- | -------- | ------------|
| miRNA-seq | 1.0, 1.1 | Quantification of known microRNAs from NGS data |
| QA        | 1.0      | Quality assessment of NGS data |
| RNA-seq   | 1.0      | Quantification of known genes from NGS data |
| WGBS      | 1.0, 1.1 | Quantification of CpG methylation from NGS data |

## Build  Docker images and push to DockerHub

The command below builds and pushes Docker images to DockerHub.
Images are tagged with tool versions.

```bash
bash build_images.sh
```

## Generate inputs files

Generating input files make WOMTOOL to validate WDL workflow and tools files.

```bash
bash generate_inputs.sh
```