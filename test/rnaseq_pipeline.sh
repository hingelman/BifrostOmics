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


######## TODO - 
###
###1-Avoid reprocessing
###Add checks to skip steps that are already done, e.g.:
###if [[ -f "$TRIMMED_DIR/${sample}_1.trimmed.fastq.gz" ]]; then
###  echo "$sample already trimmed. Skipping."
###  return
###fi
###Same for alignment and sorting.
###
###2-Better logging
###Redirect stdout/stderr to logs per step:
###
###LOG_DIR="./logs"
###mkdir -p "$LOG_DIR"
###
#### Example inside align_sample
###{
###  echo "Aligning $sample..."
###  hisat2 ...
###  samtools sort ...
###  samtools index ...
###} &> "$LOG_DIR/${sample}.align.log"
###3-add DESeq2 + plotting after featureCounts
###4-convert this to a Snakemake pipeline



# Start logging all output (stdout and stderr) to a file, *and* also print to terminal
exec > >(tee -i rnaseq_pipeline.log)
exec 2>&1

set -e  # Exit on any error
for cmd in prefetch fasterq-dump pigz fastqc multiqc fastp hisat2 samtools featureCounts parallel; do
    command -v $cmd >/dev/null 2>&1 || { echo >&2 "$cmd not found in PATH. Aborting."; exit 1; }
done

#############################
# STEP 1: Check for required input files
#############################

REFERENCE_DIR="./reference"
RAW_DIR="./raw_data"
TRIMMED_DIR="./trimmed_data"
FASTQC_DIR="./fastqc_results"
ALIGNMENTS_DIR="./alignments"
COUNTS_DIR="./counts"

mkdir -p "$RAW_DIR" "$TRIMMED_DIR" "$FASTQC_DIR" "$ALIGNMENTS_DIR" "$COUNTS_DIR"

set -eo pipefail

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

SECONDS=0
echo "Downloading SRA files..."
cd "$RAW_DIR"

SRR_IDS=$(awk -F, 'NR>1 {print $1}' runinfo.csv)

for srr in $SRR_IDS; do
  echo "Downloading $srr..."
  prefetch "$srr"
done
cd - > /dev/null || exit 1

echo "Step took $SECONDS seconds."

#############################
# STEP 3: Extract FASTQ in parallel
#############################

SECONDS=0

process_sra() {
    sra_file="$1"
    sra_dir=$(dirname "$sra_file")
    sra_base=$(basename "$sra_file" .sra)
    cd "$sra_dir" || exit 1

    # Check if both paired-end fastq.gz files exist
    if [ -f "${sra_base}_1.fastq.gz" ] && [ -f "${sra_base}_2.fastq.gz" ]; then
        echo "Skipping ${sra_base} - Paired FASTQ files already exist"
        return 0
    # Check if single-end fastq.gz exists
    elif [ -f "${sra_base}.fastq.gz" ]; then
        echo "Skipping ${sra_base} - Single-end FASTQ file already exists"
        return 0
    fi

    fasterq-dump "$(basename "$sra_file")" --threads 4

    for fq in *.fastq; do
        [ -f "$fq" ] && pigz -f -p 4 "$fq"
    done
    cd - > /dev/null || exit 1
}
export -f process_sra

echo "Extracting FASTQ files"
find "$RAW_DIR" -type f -name "*.sra" -print0 | \
    parallel -0 --jobs 4 --load 100% process_sra

echo "All .sra files processed. Step took $SECONDS seconds."

#############################
# STEP 4: FastQC + MultiQC
#############################

echo "Running FastQC and MultiQC..."
SECONDS=0

find "$RAW_DIR" -name "*.fastq.gz" | while read fq; do
    echo "Running FastQC on $fq"
 
    fastqc "$fq" -o "$FASTQC_DIR"
done

cd "$FASTQC_DIR"
multiqc .
cd - > /dev/null || exit 1

echo "FastQC and MultiQC reports created. Step took $SECONDS seconds."

#############################
# STEP 5: Trim Reads in Parallel
#############################

echo "Found $(find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l) samples for trimming"
echo "Trimming reads"
SECONDS=0

process_sample() {
    local folder="$1"
    local r1=$(find "$folder" -name "*_1.fastq.gz" | head -n1)
    local r2=$(find "$folder" -name "*_2.fastq.gz" | head -n1)
    local sample=$(basename "$folder")

    if [[ -f "$r1" && -f "$r2" ]]; then
        echo "Trimming $sample"
        echo "Input R1: $r1"
        echo "Input R2: $r2"
        echo "Output: $TRIMMED_DIR/${sample}_1.trimmed.fastq.gz"

        mkdir -p "$TRIMMED_DIR"  # to ensure it exists

        fastp \
          -i "$r1" \
          -I "$r2" \
          -o "$TRIMMED_DIR/${sample}_1.trimmed.fastq.gz" \
          -O "$TRIMMED_DIR/${sample}_2.trimmed.fastq.gz" \
          --detect_adapter_for_pe \
          --thread 2 \
          --html "$TRIMMED_DIR/${sample}.html" \
          --json "$TRIMMED_DIR/${sample}.json"
    else
        echo "Missing input files for $sample — skipping"
    fi
}
export -f process_sample
export TRIMMED_DIR
export RAW_DIR

find "$RAW_DIR" -mindepth 1 -maxdepth 1 -type d | \
    parallel --jobs $(( $(nproc) / 2 )) --load 100% process_sample

echo "All trimming jobs completed. Step took $SECONDS seconds."

#############################
# STEP 6: Genome Indexing (if needed)
#############################

if [[ ! -d "$REFERENCE_DIR/hisat2_index" ]]; then
  echo "Building HISAT2 index..."
  SECONDS=0
  mkdir -p "$REFERENCE_DIR/hisat2_index"
  hisat2-build "$FASTA" "$REFERENCE_DIR/hisat2_index/genome"
  
  echo "HISAT2 index created under /reference/hisat2_index/genome. Step took $SECONDS seconds."
  
else
  echo "HISAT2 index already exists. Skipping indexing step."
fi


#############################
# STEP 7: Align Reads to Genome
#############################

echo "Aligning reads to genome"
SECONDS=0
align_sample() {
    r1="$1"
    r2="${r1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
    sample=$(basename "$r1" | sed 's/_1\.trimmed\.fastq\.gz$//')

    # Safety check
    if [[ "$ALIGNMENTS_DIR" == "/" || -z "$ALIGNMENTS_DIR" ]]; then
        echo "ERROR: ALIGNMENTS_DIR is not set correctly!"
        exit 1
    fi

    echo "Aligning sample: $sample"
    echo "Output: $ALIGNMENTS_DIR/${sample}.sorted.bam"

    hisat2 -p 6 -x "$REFERENCE_DIR/hisat2_index/genome" \
        -1 "$r1" -2 "$r2" \
        | samtools sort -@ 2 -o "$ALIGNMENTS_DIR/${sample}.sorted.bam"

    samtools index "$ALIGNMENTS_DIR/${sample}.sorted.bam"
}

export -f align_sample
export ALIGNMENTS_DIR
export REFERENCE_DIR

find "$TRIMMED_DIR" -name "*_1.trimmed.fastq.gz" | \
        parallel --jobs $(( $(nproc) / 2 )) align_sample
        
echo "Aligned reads to genome. Step took $SECONDS seconds."

#############################
# STEP 8: Read Counting with featureCounts
#############################

if [[ ! -f "$COUNTS_DIR/gene_counts.txt" ]]; then
echo "Read Counting with featureCounts"
SECONDS=0
mkdir -p "$COUNTS_DIR"

# Find all BAM files and assign to array
mapfile -t BAM_FILES < <(find "$ALIGNMENTS_DIR" -name "*.sorted.bam")

if [ ${#BAM_FILES[@]} -eq 0 ]; then
  echo "No BAM files found in $ALIGNMENTS_DIR. Exiting."
  exit 1
fi

# Run featureCounts on all BAM files
featureCounts \
  -T 8 \
  -p \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$COUNTS_DIR/gene_counts.txt" \
  "${BAM_FILES[@]}"
  
echo "Saving sample names to $COUNTS_DIR/samples.txt"

find "$TRIMMED_DIR" -name "*_1.trimmed.fastq.gz" | sed 's|.*/||; s|_1.trimmed.fastq.gz||' > "$COUNTS_DIR/samples.txt"

echo "Read Counting with featureCounts done. Step took $SECONDS seconds"

else
  echo "featureCounts already run. Skipping."
fi

#############################

echo "RNA-seq data processing pipeline complete."