#!/bin/bash
# Location: rnaseq_project/raw_data/
# Recursively find all .sra files inside raw_data and process them
find -type f -name "*.sra" | while read sra_file; do
    echo "Processing $sra_file..."

    # Move to directory where the .sra file is, so output goes there
    sra_dir=$(dirname "$sra_file")
    cd "$sra_dir"

    # Extract fastq
    fasterq-dump "$(basename "$sra_file")" --threads 4

    # Gzip all fastq files generated in this folder
    for fq in *.fastq; do
        echo "Compressing $fq..."
        gzip -f "$fq"
    done

    # Go back to project root
    cd - > /dev/null
done

echo "✅ All .sra files processed."