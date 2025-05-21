#!/bin/bash

# Check if runinfo.csv exists
if [[ ! -f "runinfo.csv" ]]; then
  echo "runinfo.csv not found in current directory"
  exit 1
fi

# Extract SRR IDs from the Run column (skip header)
SRR_IDS=$(awk -F, 'NR>1 {print $1}' runinfo.csv)

echo "Starting download of SRR files..."

for srr in $SRR_IDS; do
  echo "Downloading $srr ..."
  prefetch "$srr"
  # Alternatively, if you want FASTQ directly:
  # fasterq-dump --split-files --gzip "$srr"
done

echo "Download finished."
