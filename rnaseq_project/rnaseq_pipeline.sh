#!/bin/bash
# File: rnaseq_pipeline.sh
# Location: rnaseq_project/ or any other name...
# Master pipeline to automate RNA-seq preprocessing and alignment

#Assumptions and readme:
# 0. Conda env bioinformaatika0 active
# 1. Structure your folders like this:
#rnaseq_project/
#├── reference/
#│   ├── genome.fna
#│   ├── annotations.gtf
#│   └── optional.gff
#├── raw_data/
#│   └── runinfo.csv
#└── rnaseq_pipeline.sh
#
#2. Place required files in reference/:
#.fna (reference genome)
#.gtf (annotation for featureCounts and alignment)
#Optional: .gff (not used by this script but can be stored)
#
#3. Create/Have runinfo.csv in raw_data/:
#Download this from SRA Run Selector
#Must contain a header and Run column with SRR IDs





set -e  # Exit on any error

#############################
# STEP 1: Check for required input files
#############################
REFERENCE_DIR="reference"
RAW_DIR="raw_data"
TRIMMED_DIR="trimmed_data"
FASTQC_DIR="fastqc_results"
ALIGNMENTS_DIR="alignments"
COUNTS_DIR="counts"

mkdir -p "$RAW_DIR" "$TRIMMED_DIR" "$FASTQC_DIR" "$ALIGNMENTS_DIR" "$COUNTS_DIR"

if [[ ! -f "$RAW_DIR/runinfo.csv" ]]; then
  echo "runinfo.csv not found in $RAW_DIR"
  exit 1
fi

FASTA=$(find "$REFERENCE_DIR" -name "*.fna" | head -n1)
GTF=$(find "$REFERENCE_DIR" -name "*.gtf" | head -n1)

if [[ ! -f "$FASTA" || ! -f "$GTF" ]]; then
  echo "Required reference files (.fna, .gtf) not found in $REFERENCE_DIR"
  exit 1
fi

#############################
# STEP 2: Download SRA files
#############################
echo "Downloading SRA files..."
cd "$RAW_DIR"
SRR_IDS=$(awk -F, 'NR>1 {print $1}' runinfo.csv)
for srr in $SRR_IDS; do
  echo "Downloading $srr..."
  prefetch "$srr"
done
cd -

#############################
# STEP 3: Extract FASTQ in parallel
#############################
process_sra() {
    sra_file="$1"
    sra_dir=$(dirname "$sra_file")
    cd "$sra_dir" || exit 1

    fasterq-dump "$(basename "$sra_file")" --threads 4

    for fq in *.fastq; do
        [ -f "$fq" ] && pigz -f -p 4 "$fq"
    done
    cd - > /dev/null || exit 1
}
export -f process_sra

find "$RAW_DIR" -type f -name "*.sra" -print0 | \
    parallel -0 --jobs 4 --load 100% process_sra

echo "All .sra files processed."

#############################
# STEP 4: FastQC + MultiQC
#############################
find "$RAW_DIR" -name "*.fastq.gz" | while read fq; do
    echo "Running FastQC on $fq"
    fastqc "$fq" -o "$FASTQC_DIR"
done

cd "$FASTQC_DIR"
multiqc .
cd -

#############################
# STEP 5: Trim Reads in Parallel
#############################
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
          -o "$TRIMMED_DIR/${sample}_1.trimmed.fastq.gz" \
          -O "$TRIMMED_DIR/${sample}_2.trimmed.fastq.gz" \
          --detect_adapter_for_pe \
          --thread 2 \
          --html "$TRIMMED_DIR/${sample}.html" \
          --json "$TRIMMED_DIR/${sample}.json"
    fi
}
export -f process_sample

find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type d | \
    parallel --jobs 4 --load 100% process_sample

echo "All trimming jobs completed."

#############################
# STEP 6: Genome Indexing (if needed)
#############################
if [[ ! -d "$REFERENCE_DIR/hisat2_index" ]]; then
  echo "Building HISAT2 index..."
  mkdir -p "$REFERENCE_DIR/hisat2_index"
  hisat2-build -p 6 "$FASTA" "$REFERENCE_DIR/hisat2_index/genome"
fi

#############################
# STEP 7: Align Reads to Genome
#############################
align_sample() {
    r1="$1"
    r2="${r1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
    sample=$(basename "$r1" | cut -d_ -f1)

    hisat2 -p 6 -x "$REFERENCE_DIR/hisat2_index/genome" \
        -1 "$r1" -2 "$r2" \
        | samtools sort -@ 2 -o "$ALIGNMENTS_DIR/${sample}.sorted.bam"

    samtools index "$ALIGNMENTS_DIR/${sample}.sorted.bam"
}
export -f align_sample

find "$TRIMMED_DIR" -name "*_1.trimmed.fastq.gz" | \
    parallel --jobs 1 align_sample

#############################
# STEP 8: Read Counting with featureCounts
#############################
BAM_FILES=$(find "$ALIGNMENTS_DIR" -name "*.sorted.bam" | tr '\n' ' ')
featureCounts -p -T 8 -a "$GTF" -o "$COUNTS_DIR/gene_counts.txt" $BAM_FILES

echo "RNA-seq pipeline complete."
