#!/bin/bash

#SBATCH --partition=workq
#SBATCH --job-name=01_NextITS
#SBATCH --output=./01_NextITS.out
#SBATCH --error=./01_NextITS.err
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --mem=36G


module load bioinfo/Nextflow/25.04.0
module load containers/singularity/3.9.9 
module load devel/java/17.0.6

#To get the latest version of NextITS from the GitHub 
nextflow pull vmikk/NextITS

#To get the nextits.sif file to load all the module and packages to run NextITS
singularity pull nextits.sif docker://vmikk/nextits:1.1.0

#To get the help message with all the option and their meanings
nextflow run vmikk/NextITS --help

#Change the primer if it is not the one used for the full ITS region and also adapt the --its_region
#Not mandatory to use a chimera database, only denovo can be done.
nextflow run vmikk/NextITS -r main \
  -resume \
  --demultiplexed true \
  --step "Step1" \
  --input "$(pwd)/path/to/fastq_files" \
  --primer_forward GTACACACCGCCCGTCG \
  --primer_reverse CGCCTSCSCTTANTDATATGC \
  --its_region "full" \
  --hp "true" \
  --chimera_methods "ref,denovo" \
  --chimera_db "$(pwd)/path/to/chimera_database" \
  --tj "true" \
  --outdir  "Step1_Results/test" \
  -with-singularity "$(pwd)/nextits.sif"
