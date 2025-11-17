#!/bin/bash
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=30G

# Requirements:
# - vsearch >= 2.29
# - mumu (https://github.com/frederic-mahe/mumu.git)
# - gcc >= 10.
#   e.g -> install via: git clone https://github.com/frederic-mahe/mumu.git && cd mumu && module load compilers/gcc/12.2.0 && make
# - To check: in the mumu directory do, ./mumu --help

#Computation parameters:
 #--minimum_match FLOAT                minimum similarity threshold (84.0)
 #--minimum_ratio FLOAT                minimum abundance ratio (1.0)
 #--minimum_ratio_type STRING          "min" or "avg" abundance ratio ("min")
 #--minimum_relative_cooccurence FLOAT relative father-son spread (0.95)

#module purge
module load bioinfo/VSEARCH/2.29.3
module load compilers/gcc/12.2.0



IDENTITY_PERCENTAGE=0.95
OUTPUT_DIR=LULU_OTU_CLUSTERING

while getopts "f:s:o:t:c:" option
do
        case $option in
                h)
                    echo "Usage: $0 -f <fasta_file> -o <output_directory> -t <otu_table> -s <simmilarity_cutoff>"
                    echo "  -f  fasta file (can be compressed)"
                    echo "  -o  Output directory (default: LULU_OTU_CLUSTERING)"
                    echo "  -s  Identity percentage for curation (default: 0.95)"
                    echo "  -h  Show this help message"
		    echo "  -t  OTU table (OTU_name tab_sep sample1... ) (can be compressed)"
		    echo "  -c  OTU table complete (OTU_name abundance length... sample1... ) from the last step (can be compressed)"
                    exit 0
                    ;;
                f)
                    FASTA_FILE="$OPTARG"
                    ;;
                o)
                    OUTPUT_DIR="$OPTARG"
                    ;;
                s)
                    IDENTITY_PERCENTAGE="$OPTARG"
                    ;;
                t)
                    OTU_TABLE="$OPTARG"
                    ;;
		c)
		    COMPLETE_OTU_TABLE="$OPTARG"
		    ;;
        esac
done

DIR="$(pwd)"
cd $DIR
mkdir -p $OUTPUT_DIR

# First step: Build the match list file
vsearch \
      --usearch_global ${FASTA_FILE} \
      --db             ${FASTA_FILE} \
      --self \
      --id         ${IDENTITY_PERCENTAGE} \
      --iddef      1 \
      --gapopen 20I/2E \
      --gapext  2I/1E \
      --query_cov  0.9 \
      --userfields query+target+id \
      --maxaccepts 0 \
      --maxhits    0 \
      --threads    8 \
      --userout    $OUTPUT_DIR/LULU_match_list.txt

#To get only the name and not the name;size=... in ordre to match the name on your OTU table
sed -E 's/;size=[0-9]+//g; s/;//g' $OUTPUT_DIR/LULU_match_list.txt > $OUTPUT_DIR/LULU_match_list_clean.txt
mv $OUTPUT_DIR/LULU_match_list_clean.txt $OUTPUT_DIR/LULU_match_list.txt

#Second step: Running MUMU

#For compressed OTU table (*.gz)
if [[ "$OTU_TABLE" == *.gz ]]; then
	gunzip --stdout "$OTU_TABLE" > "$OUTPUT_DIR/tmp_OTU_table.txt"
else
	cp $OTU_TABLE  $OUTPUT_DIR/tmp_OTU_table.txt
fi


MATCH=$(awk -v a=${IDENTITY_PERCENTAGE} 'BEGIN {print(a*100) }')

#Same parameters used in NextITS pipeline
mumu/mumu \
      --otu_table     $OUTPUT_DIR/tmp_OTU_table.txt \
      --match_list    $OUTPUT_DIR/LULU_match_list.txt \
      --new_otu_table $OUTPUT_DIR/OTU_table_LULU.txt \
      --log 	      $OUTPUT_DIR/LULU_merging_statistics.txt \
      --threads                      8 \
      --minimum_match                ${MATCH} \
      --minimum_ratio                1 \
      --minimum_ratio_type           "min" \
      --minimum_relative_cooccurence 0.95

#To get the full OTU_TABLE after LULU

if [[ "$COMPLETE_OTU_TABLE" == *.gz ]]; then
        gunzip --stdout "$COMPLETE_OTU_TABLE" > "$OUTPUT_DIR/tmp_OTU_table_full.txt"
else
        cp $COMPLETE_OTU_TABLE  $OUTPUT_DIR/tmp_OTU_table_full.txt
fi
awk 'NR==FNR{
    # lire OTU_table_OTU_filtered.txt pour stocker infos par amplicon
    if(FNR>1){
        key=$3        # colonne amplicon
        abundance[key]=$2
        lengtha[key]=$4
        chimera[key]=$5
        spread[key]=$6
        identity[key]=$7
        taxonomy[key]=$8
    }
    next
}
FNR==1{
    # entête du fichier OTU_table_LULU.txt
    print "amplicon\tabundance\tlength\tchimera\tspread\tidentity\ttaxonomy\t"$0
    next
}
{
    key=$1
    # rajouter infos du OTU_filtered si existantes
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", key, abundance[key], lengtha[key], chimera[key], spread[key], identity[key], taxonomy[key]
    # rajouter colonnes samples du fichier LULU
    for(i=2;i<=NF;i++) printf "\t%s",$i
    print ""
}'  $OUTPUT_DIR/tmp_OTU_table_full.txt $OUTPUT_DIR/OTU_table_LULU.txt > $OUTPUT_DIR/OTU_table_LULU_full.txt

rm $OUTPUT_DIR/tmp*
