#!/bin/bash
set -e  # Exit on error

######################

# Prerequisites:
# 1. Conda installed (https://docs.conda.io/)
# 2. Your project directory structure is like this:
# folder/
#├── reference/			---- reference genome info here
#│   ├── genome.fna
#│   └── annotations.gtf
#├── raw_data/
#│   └── runinfo.csv     ---- this file specifies the samples
#├── Snakefile
#├── config/
#│   ├── config.yaml
#├── envs/
#│   └── rnaseq.yaml
#└── run_Snakemake_workflow.sh

######################

# 1. Creating (if needed) and activating conda environment
echo "=== Creating or activating conda environment ==="
if conda env list | grep -q "rnaseq_pipeline"; then
    echo "Environment exists, skipping creation."
else
    conda env create -f envs/rnaseq.yaml -n rnaseq_pipeline
fi

# enabling 'conda activate' in shell scripts
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnaseq_pipeline

# 2. Test individual tools
echo "=== Testing tool installations ==="
tools=(
    "snakemake --version"
    "prefetch --version"
    "fasterq-dump --version"
    "fastqc --version"
    "multiqc --version"
    "fastp --version"
    "hisat2 --version"
    "samtools --version"
    "featureCounts -v"
)

for cmd in "${tools[@]}"; do
    echo "Testing: $cmd"
    $cmd || echo "Warning: $cmd failed"
done

# 3. Verify config and input files
echo "=== Verifying config and input files ==="
test -f config/config.yaml || { echo "Error: config.yaml missing"; exit 1; }
test -f raw_data/runinfo.csv || { echo "Error: runinfo.csv missing"; exit 1; }

# 4. Run the pipeline
echo "=== Running the Snakemake workflow ==="
snakemake --use-conda --cores 8 --printshellcmds --conda-frontend conda --rerun-incomplete --latency-wait 10