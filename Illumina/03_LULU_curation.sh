#!/bin/bash

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


NB_CORES=1

while getopts "hf:o:s:t:c:m:n:" option
do
        case $option in
                h)
                    echo "Usage: $0 -f <fasta_file> -o <output_directory> -t <otu_table> -s <similarity_cutoff>"
                    echo "  -f  Fasta file (can be compressed)"
                    echo "  -o  Output directory (default: 03_LULU_OTU_CLUSTERING)"
                    echo "  -s  Identity percentage for curation (default: 0.95)"
                    echo "  -h  Show this help message"
                    echo "  -m  Command to run MUMU (default: mumu)"
		    		echo "  -t  OTU table without metadata(OTU_name tab_sep sample1... ) (can be compressed)"
		   	 		echo "  -c  OTU table complete with metadata (OTU_name abundance length... sample1... ) from the last step (can be compressed)"
                    echo "  -n  Number of CPU cores to use (default = 1)"
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
             	m)
                    MUMU="$OPTARG"
                    ;;
				c)
		    		COMPLETE_OTU_TABLE="$OPTARG"
		    		;;
		    	n)
		    		NB_CORES="$OPTARG"
		    		;;
        esac
done


if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$(pwd)/03_LULU_OTU_CLUSTERING"
    mkdir -p "$OUTPUT_DIR"
fi

if [ -z "$MUMU" ]; then
    MUMU="mumu"
fi

if [ -z "$IDENTITY_PERCENTAGE" ]; then
    IDENTITY_PERCENTAGE=0.95
fi


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
      --threads    $NB_CORES \
      --userout    $OUTPUT_DIR/LULU_match_list.txt

# To get only the name and not the name;size=... in order to match the name on the OTU table
sed -E 's/;size=[0-9]+//g; s/;//g' $OUTPUT_DIR/LULU_match_list.txt > $OUTPUT_DIR/LULU_match_list_clean.txt
mv $OUTPUT_DIR/LULU_match_list_clean.txt $OUTPUT_DIR/LULU_match_list.txt


# Second step: Running MUMU

# For compressed OTU table (*.gz)
if [[ "$OTU_TABLE" == *.gz ]]; then
	gunzip --stdout "$OTU_TABLE" > "$OUTPUT_DIR/tmp_OTU_table.txt"
else
	cp $OTU_TABLE  $OUTPUT_DIR/tmp_OTU_table.txt
fi


MATCH=$(awk -v a=${IDENTITY_PERCENTAGE} 'BEGIN {print(a*100) }')

$MUMU \
      --otu_table     $OUTPUT_DIR/tmp_OTU_table.txt \
      --match_list    $OUTPUT_DIR/LULU_match_list.txt \
      --new_otu_table $OUTPUT_DIR/OTU_table_LULU.txt \
      --log 	      $OUTPUT_DIR/LULU_merging_statistics.txt \
      --minimum_match                ${MATCH} \
      --minimum_ratio                1 \
      --minimum_ratio_type           "min" \
      --minimum_relative_cooccurence 0.95

# To get the full OTU_TABLE after LULU:

if [[ "$COMPLETE_OTU_TABLE" == *.gz ]]; then
        gunzip --stdout "$COMPLETE_OTU_TABLE" > "$OUTPUT_DIR/tmp_OTU_table_full.txt"
else
        cp $COMPLETE_OTU_TABLE  $OUTPUT_DIR/tmp_OTU_table_full.txt
fi
awk 'NR==FNR{
    # red OTU_table_OTU_filtered.txt to stock infos per amplicon
    if(FNR>1){
        key=$3        # amplicon column
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
    # header of file OTU_table_LULU.txt
    print "amplicon\tabundance\tlength\tchimera\tspread\tidentity\ttaxonomy\t"$0
    next
}
{
    key=$1
    # add infos of OTU_filtered if exist
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s", key, abundance[key], lengtha[key], chimera[key], spread[key], identity[key], taxonomy[key]
    # add columns samples of LULU file
    for(i=2;i<=NF;i++) printf "\t%s",$i
    print ""
}'  $OUTPUT_DIR/tmp_OTU_table_full.txt $OUTPUT_DIR/OTU_table_LULU.txt > $OUTPUT_DIR/tmp2_OTU_table_LULU_full.txt

head -1 $OUTPUT_DIR/tmp2_OTU_table_LULU_full.txt > $OUTPUT_DIR/tmp_header
sed 1d  $OUTPUT_DIR/tmp2_OTU_table_LULU_full.txt | sort -k2 -rn > $OUTPUT_DIR/tmp3_OTU_table_LULU_full.txt

cat $OUTPUT_DIR/tmp_header $OUTPUT_DIR/tmp3_OTU_table_LULU_full.txt > $OUTPUT_DIR/OTU_table_LULU_full.txt

rm $OUTPUT_DIR/tmp*
