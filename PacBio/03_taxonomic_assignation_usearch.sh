#!/bin/bash

#SBATCH --partition=workq
#SBATCH --job-name=assignation
#SBATCH --output=./logs/assignation.out
#SBATCH --error=./logs/assignation.err
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --mem=12G

module load bioinfo/VSEARCH/2.29.3

#Assignation wth usearch_global
vsearch --usearch_global Path/to/fasta_file_after_step_2 \
  	--db Path/to/the/database  \
   	--id 0.7    \
	--iddef 2 \
   	--strand plus  \
   	--threads 8 \
	--dbmask none \
        --qmask none \
        --rowlen 0 \
        --notrunclabels \
        --userfields query+id+target \
        --maxaccepts 1 \
        --maxrejects 32 \
        --userout path/to/the/output_directory
