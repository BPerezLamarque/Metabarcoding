#!/bin/bash

#SBATCH --partition=workq
#SBATCH --job-name=02_NextITS
#SBATCH --output=./logs/02_NextITS.out
#SBATCH --error=./logs/02_NextITS.err
#SBATCH --cpus-per-task=10
#SBATCH --time=96:00:00
#SBATCH --mem=36G


module load bioinfo/Nextflow/25.04.0
module load containers/singularity/3.9.9 
module load devel/java/17.0.6


nextflow run vmikk/NextITS -r main \
  --step "Step2" \
  --data_path  "$(pwd)/Step1_Results/" \
  --preclustering "unoise" \
  --clustering "vsearch" \
  --otu_id 0.97 \
  --otu_iddef 2 \
  --lulu "true" \
  --outdir     "Step2_vsearch" \
  -with-singularity "$(pwd)/nextits.sif"
