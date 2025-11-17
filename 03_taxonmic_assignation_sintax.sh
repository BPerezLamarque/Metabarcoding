#!/bin/bash

#SBATCH --partition=workq
#SBATCH --job-name=assignation
#SBATCH --output=./logs/sintax_assignation.out
#SBATCH --error=./logs/sintax_assignation.err
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --mem=12G

module load bioinfo/VSEARCH/2.29.3

#Sintax assignation
vsearch --sintax Path/to/the/fasta_file_from_step_2 \
  	--db Path/to/the/database  \
   	--sintax_cutoff 0.8    \
   	--tabbedout Path/to/the/output_directory_file \
   	--strand plus    \
   	--threads 8

t1=$(date "+%s")

echo "elapsed wall time is $((t1-t0))"
