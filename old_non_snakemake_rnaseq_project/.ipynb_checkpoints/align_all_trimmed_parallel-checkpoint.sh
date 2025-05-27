#!/bin/bash
# File: align_all_trimmed_parallel.sh
# Purpose: Align all paired-end trimmed FASTQ files in parallel using HISAT2 + samtools
# Hardware: 4 cores, 8 threads total available

mkdir -p bam_files

# Function to align one sample
align_sample() {
    local sample="$1"
    local f1="trimmed_data/${sample}_1.trimmed.fastq.gz"
    local f2="trimmed_data/${sample}_2.trimmed.fastq.gz"
    local sam_out="bam_files/${sample}.sam"
    local bam_out="bam_files/${sample}.bam"

    echo "Aligning $sample ..."

    # Use 6 threads for HISAT2 (leaving 2 threads for other processes)
    hisat2 -x reference/ca_index \
           -1 "$f1" \
           -2 "$f2" \
           -p 6 \
           -S "$sam_out"

    # Use 2 threads for samtools conversion
    samtools view -@ 2 -S -b "$sam_out" > "$bam_out"
    rm "$sam_out"
}

export -f align_sample

# Get sample names and run in parallel
# Run only 1 alignment at a time since HISAT2 will use 6 threads
find trimmed_data -name "*_1.trimmed.fastq.gz" | sed 's/.*\///; s/_1.trimmed.fastq.gz//' | \
    parallel --jobs 1 --load 100% align_sample

echo "All alignments complete!"