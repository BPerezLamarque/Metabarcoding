#!/bin/bash
#SBATCH --job-name=trim_db
#SBATCH --output=trim_db.out
#SBATCH --error=trim_db.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module purge
module load bioinfo/Cutadapt/5.0

# Paths
WORKDIR="$HOME/work/Stage_2A/donnees_miti/20250919_DEFIS/hifi_reads"
INPUT="${WORKDIR}/database/SINTAX_EUK_ITS_v1.9.4.fasta"
OUTPUT_DIR="$HOME/work/Stage_2A/donnees_illumina"

mkdir -p "$OUTPUT_DIR"

bash cut_database.sh \
    -i "$INPUT" \
    -f GGAAGTAAAAGTCGTAACAAGG \
    -r CAAGAGATCCGTTGTTGAAAGTK \
    -o "$OUTPUT_DIR" \
    -n General_EUK_ITS_v1.9.4_Ted \
    -x false \
    -y true \
    -l 50 \
    -t "sintax"

