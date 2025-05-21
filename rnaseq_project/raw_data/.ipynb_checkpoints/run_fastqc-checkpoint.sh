#!/bin/bash
# File: run_fastqc.sh
# Location: rnaseq_project/raw_data/

mkdir -p ../fastqc_results

find . -name "*.fastq.gz" | while read fq; do
    echo "Running FastQC on $fq"
    fastqc "$fq" -o ../fastqc_results
done