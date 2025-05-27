#!/bin/bash
# File: trim_all_fastq_parallel.sh
# Location: rnaseq_project/raw_data/
# Optimized for 8 threads across 4 cores

# Create output directory
mkdir -p ../trimmed_data

# Function to process a single sample
process_sample() {
    local folder="$1"
    local r1=$(find "$folder" -name "*_1.fastq.gz" | head -n1)
    local r2=$(find "$folder" -name "*_2.fastq.gz" | head -n1)
    local sample=$(basename "$folder")

    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "Trimming $sample"
        
        fastp \
          -i "$r1" \
          -I "$r2" \
          -o ../trimmed_data/"${sample}"_1.trimmed.fastq.gz \
          -O ../trimmed_data/"${sample}"_2.trimmed.fastq.gz \
          --detect_adapter_for_pe \
          --thread 2 \
          --html ../trimmed_data/"${sample}".html \
          --json ../trimmed_data/"${sample}".json
    fi
}

export -f process_sample

# Find all directories and process them in parallel
find . -maxdepth 1 -type d ! -path . | \
    parallel --jobs 4 --load 100% process_sample

echo "All trimming jobs completed!"