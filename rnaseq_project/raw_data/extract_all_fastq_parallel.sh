#!/bin/bash
# Location: rnaseq_project/raw_data/
# File: extract_all_fastq_parallel.sh
# Recursively find all .sra files inside raw_data and process them
# 
# 
# First, install required dependencies:
# On Ubuntu/Debian
#   sudo apt-get install parallel
#   sudo apt-get install parallel pigz
# 
# 
# Function to process a single SRA file
process_sra() {
    sra_file="$1"
    echo "Processing $sra_file..."

    # Move to directory where the .sra file is
    sra_dir=$(dirname "$sra_file")
    cd "$sra_dir" || exit 1

    # Extract fastq using fasterq-dump
    # Using more threads per file since we have 8 total threads
    fasterq-dump "$(basename "$sra_file")" --threads 4

    # Gzip all fastq files generated in this folder using pigz
    for fq in *.fastq; do
        if [ -f "$fq" ]; then
            echo "Compressing $fq..."
            pigz -f -p 4 "$fq"  # Using 4 threads for compression
        fi
    done

    cd - > /dev/null || exit 1
}

# Export the function so GNU parallel can use it
export -f process_sra

# Find all .sra files and process them in parallel
# Process more files simultaneously to maximize thread usage
find . -type f -name "*.sra" -print0 | \
    parallel -0 --jobs 4 --load 100% process_sra

echo "✅ All .sra files processed."