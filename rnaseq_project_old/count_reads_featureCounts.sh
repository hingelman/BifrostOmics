#!/bin/bash
# File: count_reads_featureCounts.sh
# Location: rnaseq_project folder
# Description: Count gene-level reads from BAMs

mkdir -p counts

featureCounts \
  -T 8 \
  -p \
  -t gene \
  -g gene_id \
  -a reference/GCF_000484505.2_ASM48450v2_genomic.gtf \
  -o counts/gene_counts.txt \
  bam_files/*.bam