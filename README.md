BifrostOmics ver 0.1 This project will aim to be a bridge between omics datasets, bioinformatics & non-CS users by providing a collection of reproducible pipelines for various omics datasets:

RNA-seq Pipeline
placeholder1
placeholder2
placeholder3
etc...


#
pipeline specific readmes:

#
#
# # RNA-seq Pipeline v0.1: Differential Gene Expression from your dataset
This workflow implements a reproducible RNA-seq pipeline using publicly available tools. It includes quality control, trimming, alignment, read summarization, and automation via Snakemake.

#Pipeline includes:

Download SRA/FASTQ files based on runinfo.csv
QC with FastQC and MultiQC
Trim reads with fastp
Build HISAT2 index and align reads with hisat2
Index/Align BAM files
Count reads with featureCounts
Prerequisites:
1. Conda installed (https://docs.conda.io/)
2. Your project directory structure is like this:
folder/
#├── reference/ ---- reference genome info here 
#│ ├── genome.fna 
#│ └── annotations.gtf 
#├── raw_data/ 
#│ └── runinfo.csv ---- this file specifies the samples 
#├── Snakefile 
#├── config/ 
#│ ├── config.yaml 
#├── envs/ 
#│ └── rnaseq.yaml 
#└── run_Snakemake_workflow.sh

Running
Clone the folder "rnaseq_project"
Run the script "run_Snakemake_workflow.sh" All software is installed via Conda using envs/*.yaml which is embedded in "run_Snakemake_workflow.sh" script To visualize the DAG add the "--dag", eg "bash run_Snakemake_workflow.sh --dag"
Outputs
fastqc_results/: Quality reports
trimmed_data/: Trimmed FASTQ files
alignments/: BAM files
counts/: Gene counts file(s)
Tools Used
conda
snakemake
prefetch
fasterq-dump
fastqc
multiqc
fastp
hisat2
samtools
featureCounts
#