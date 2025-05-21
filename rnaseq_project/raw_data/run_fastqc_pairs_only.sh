#!/bin/bash

# Exit on error
set -e

# Find all _1.fastq or _1.fastq.gz files
find . -type f \( -name "*_1.fastq" -o -name "*_1.fastq.gz" \) | sort | while read r1; do
    # Derive the R2 filename by replacing _1 with _2
    r2="${r1/_1./_2.}"

    # Check if R2 file exists
    if [[ -f "$r2" ]]; then
        echo "Running FastQC on: $r1 and $r2"
        fastqc "$r1" "$r2"
    else
        echo "Warning: Paired file not found for $r1"
    fi
done

