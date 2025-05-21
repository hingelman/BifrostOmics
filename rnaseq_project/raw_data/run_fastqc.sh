#!/bin/bash

# Exit on error
set -e

# Find all fastq or fastq.gz files in subdirectories
find . -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort | while read file; do
    echo "Running FastQC on $file"
    fastqc "$file"
done