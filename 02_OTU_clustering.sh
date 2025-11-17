#!/bin/bash

#SBATCH --partition=workq
#SBATCH --output=./02.out
#SBATCH --error=./02.err
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --mem=36G

# This script performs OTU clustering using VSEARCH and generates a contingency table.
# It requires the following software: VSEARCH 2.29.3, Python 3.11.1.
# It also requires the following parameters: path to the directory containing input fasta files, path to the database, path to the scripts, output directory, abundance threshold, and identity percentage.

module purge
module load bioinfo/VSEARCH/2.29.3
module load devel/python/Python-3.11.1
module load bioinfo/swarm/3.1.3
module load bioinfo/SeqKit/2.9.0
module load statistics/R/4.4.0

nb_cores=2
OUTPUT_DIR=OTU_CLUSTERING
IDENTITY_PERCENTAGE=0.97
CUTOFF_PROB=0.01
path_scripts=$(pwd)/python_scripts
METHOD="usearch"
RESUME="true"
UNOISE="true"
CLUSTER="vsearch"
TAGJUMP="true"

while getopts "i:s:d:o:c:m:p:r:u:v:t:" option
do
        case $option in
                h)
                    echo "Usage: $0 -i <input_files> -s <path_scripts> -d <path_db> -o <output_directory> -c <identity_percentage> -m <method> -p <proba_cutoff> -r <resume> -v <cluster_method> -t <tagjump>"
		    echo "E.g: sbatch 02_OTU_clustering.sh -i Path_to/Demultiplexed_data/ -d Path_to_database/*.fasta -c 0.97 -o OTU_CLUSTERING -p 1 -u "false" -m "sintax" "
                    echo "  -i  Input directory containing fasta files (can be compressed)"
                    echo "  -d  Path to the database"
                    echo "  -o  Output directory (default: OTU_CLUSTERING)"
                    echo "  -c  Identity percentage for clustering (default: 0.97)"
                    echo "  -h  Show this help message"
		    echo "  -m  method for taxonomic assignation: sintax or usearch (default: "usearch")"
		    echo "  -p  Cutoff for sintax probability for assignation (default: 0.0)"
		    echo "  -r  resume where the last outpute file was created (default: true)"
		    echo "  -u  Do a denoising step with unoise (default = true)"
		    echo "  -s  Path to your python scripts"
		    echo "  -v  Clustering method: vsearch or swarm (default = "vsearch")"
		    echo "  -t  Do a TagJump filtration (default = true)"
                    exit 0
                    ;;
                i)
		    INPUT_FILES="$OPTARG"
                    ;;
                d)
                    path_to_db="$OPTARG"
                    ;;
                o)
                    OUTPUT_DIR="$OPTARG"
                    ;;    
                c)
                    IDENTITY_PERCENTAGE="$OPTARG"
                    ;;
		m)
		    METHOD="$OPTARG"
	       	    ;;
		p)
		    CUTOFF_PROB="$OPTARG"
		    ;;
		r)
		    RESUME="$OPTARG"
		    ;;
		u)
		    UNOISE="$OPTARG"
		    ;;
		v)
		    CLUSTER="$OPTARG"
		    ;;
		t)
		    TAGJUMP="$OPTARG"
		    ;;
        esac
done

if [ -v "$path_to_db" ]; then
    echo "Warning: No database provided. Skipping taxonomic assignment."
fi

DIR="$(pwd)"
cd $DIR
mkdir -p $OUTPUT_DIR

#STEP 2
#OTU Clustering

# indicate the name of the database for taxonomic assignation (generated in Step 0):

if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_amplicon_derep.fasta" ] && [ -s "$OUTPUT_DIR/reads_amplicon_derep.uc" ]; then
        echo "Output file already exists, skipping this step."
else
	echo "Running Dereplication by vsearch on: $OUTPUT_DIR/reads_amplicon.fasta"
	cat $INPUT_FILES/*.f* > $OUTPUT_DIR/reads_amplicon.fasta
	vsearch \
		--derep_fulllength $OUTPUT_DIR/reads_amplicon.fasta \
    		--sizein \
    		--sizeout \
    		--relabel_sha1 \
    		--fasta_width 0 \
		--uc $OUTPUT_DIR/reads_amplicon_derep.uc \
    		--output $OUTPUT_DIR/reads_amplicon_derep.fasta
fi

###TAG-JUMP REMOVAL
#NEED dereplication.parquet and SeqTab.txt with sampleID seqID abundance

if [ "$TAGJUMP" = "true" ]; then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/Seq_tab_TagJumpFiltered.fasta" ]; then
                echo "Output file already exists, skipping this step."
	else
		for f in "$INPUT_FILES"/*.f*; do
		SAMPLE=$(basename "$f" | sed -E 's/(\.[^.]+)+$//')
		cat "$f" \
			| sed -r '/^>/ s/;sample=[^;]*//g; s/;;/;/g' \
			| sed "s/^>.*/&;sample=$SAMPLE;/ ; s/_NoChimera//g ; s/_RescuedChimera//g ; s/_JoinedPE//g ; s/_Homopolymer_compressed//g"
		done \
			| vsearch --sortbysize - --sizein --sizeout --fasta_width 0 --output - \
			| sed -r '/^>/ s/;;/;/g' >  $OUTPUT_DIR/Seq_not_filtered.fasta

		seqkit seq --name $OUTPUT_DIR/Seq_not_filtered.fasta \
			| sed 's/;/\t/g; s/size=//; s/sample=// ; s/\t*$//' \
			| awk -F"\t" '{print $2 "\t" $1 "\t" $3}' > $OUTPUT_DIR/Seq_tab_not_filtered.txt

		awk '$1=="S" {print $9 "\t" $9} $1=="H" {print $9 "\t" $10}' $OUTPUT_DIR/reads_amplicon_derep.uc \
			> $OUTPUT_DIR/reads_amplicon_derep.tsv
		sed -E 's/;size=[0-9]+//g; s/;//g' $OUTPUT_DIR/reads_amplicon_derep.tsv |sort > $OUTPUT_DIR/reads_amplicon_derep_filtered.tsv
		mv $OUTPUT_DIR/reads_amplicon_derep_filtered.tsv $OUTPUT_DIR/reads_amplicon_derep.tsv

		Rscript $path_scripts/tag_jump_removal_longtab.R \
    			--seqtab $OUTPUT_DIR/Seq_tab_not_filtered.txt \
    			--uc $OUTPUT_DIR/reads_amplicon_derep.tsv \
    			--uncross_f 0.01 \
    			--uncross_p 1 \
    			--outdir $OUTPUT_DIR

		awk '{if(NR>1) print ">"$2";sample="$1"\t"$3}' $OUTPUT_DIR/Seq_tab_TagJumpFiltered.txt > $OUTPUT_DIR/SeqID_sample_to_abund.txt

		awk 'NR==FNR{map=$1;abundance[map]=$2;next} /^>/{split($1,s,";size="); key=s[1]; ab=s[2]; if(key in abundance){print key";size="abundance[key]; keep=1} else {keep=0}; next} {if(keep) print}' $OUTPUT_DIR/SeqID_sample_to_abund.txt $OUTPUT_DIR/Seq_not_filtered.fasta > $OUTPUT_DIR/Seq_tab_TagJumpFiltered.fasta

		vsearch \
  			--derep_fulllength $OUTPUT_DIR/Seq_tab_TagJumpFiltered.fasta \
  			--sizein \
  			--sizeout \
  			--fasta_width 0 \
  			--relabel_sha1 \
  			--output $OUTPUT_DIR/Seq_tab_TagJumpFiltered_derep.fasta \
  			--uc $OUTPUT_DIR/Seq_tab_TagJumpFiltered_derep.uc
		mv $OUTPUT_DIR/Seq_TagJumpFiltered_derep.fasta $OUTPUT_DIR/reads_amplicon_derep.fasta
	fi
else
	 echo "Skipping the TagJump step."
fi


if [ "$UNOISE" = "true" ]; then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/UNOISE.uc" ]; then
		echo "Output file already exists, skipping this step."
	else
		echo "Running UNOISE by vsearch on: $OUTPUT_DIR/reads_amplicon_derep.fasta"
		vsearch \
			--cluster_unoise $OUTPUT_DIR/reads_amplicon_derep.fasta \
      			--unoise_alpha  6  \
      			--minsize 1 \
      			--iddef   2 \
      			--qmask   dust \
      			--gapopen  20I/2E \
      			--gapext   2I/1E \
      			--threads 8 \
      			--fasta_width 0 \
      			--sizein --sizeout \
      			--centroids $OUTPUT_DIR/reads_amplicon_UNOISE.fa \
      			--uc $OUTPUT_DIR/UNOISE.uc
	mv $OUTPUT_DIR/reads_amplicon_UNOISE.fa $OUTPUT_DIR/reads_amplicon_derep.fasta
	fi
else
	echo "Skipping the unoise step."
fi


if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_amplicon_sorted.fasta" ]; then
	echo "Output file already exists, skipping this step."
else
	echo "Running sortbysize by vsearch on: $OUTPUT_DIR/reads_amplicon_derep.fasta"
	vsearch -sortbysize $OUTPUT_DIR/reads_amplicon_derep.fasta -output $OUTPUT_DIR/reads_amplicon_sorted.fasta -minsize 1
	
	rm $OUTPUT_DIR/reads_amplicon_derep.fasta
fi

if [ "$CLUSTER" = "vsearch" ];then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_OTU97_final.fasta" ] &&  [ -s "$OUTPUT_DIR/clusters_OTU97.uc" ] && [ -s "$OUTPUT_DIR/reads_mapped_OTU97.txt" ] &&  [ -s "$OUTPUT_DIR/stats_file_OTU97.txt" ]; then
		echo "Output file already exists, skipping this step."
	else
		echo "Running vsearch clusterisation on: $OUTPUT_DIR/reads_amplicon_sorted.fasta"
		vsearch -cluster_size  $OUTPUT_DIR/reads_amplicon_sorted.fasta \
			--threads $nb_cores \
    			--id $IDENTITY_PERCENTAGE --centroids $OUTPUT_DIR/reads_OTU97.fasta \
    			--uc $OUTPUT_DIR/clusters_OTU97.uc \
    			--sizein --sizeout
	
		vsearch --fasta_width 0 --sortbysize $OUTPUT_DIR/reads_OTU97.fasta --output $OUTPUT_DIR/reads_OTU97_final.fasta

		python3 $path_scripts/map2qiime.py $OUTPUT_DIR/clusters_OTU97.uc > $OUTPUT_DIR/reads_mapped_OTU97.txt
		python3 $path_scripts/make_stats.py $OUTPUT_DIR/reads_OTU97_final.fasta > $OUTPUT_DIR/stats_file_OTU97.txt
	fi
elif [ "$CLUSTER" = "swarm" ];then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/SWARM_representatives.fasta" ] &&  [ -s "$OUTPUT_DIR/SWARM.uc" ] && [ -s "$OUTPUT_DIR/SWARM.struct" ] &&  [ -s "$OUTPUT_DIR/SWARM.stats" ]; then
                echo "Output file already exists, skipping this step."
	else
		echo "Running SWARM clusterisation on: $OUTPUT_DIR/reads_amplicon_sorted.fasta"
		
		awk '/^>/ {if (seq) print header "\n" seq; header=$0; seq=""} 
		!/^>/ {gsub(/[ \t\r\n]/,""); seq=seq $0} 
		END {if (seq) print header "\n" seq}' $OUTPUT_DIR/reads_amplicon_sorted.fasta > $OUTPUT_DIR/reads_amplicon_uniline.fasta
		awk '/^>/ {if (seq && seq !~ /[^ACGT]/) print header "\n" seq; header=$0; seq=""}
     		!/^>/ {seq=$0}
     		END {if (seq && seq !~ /[^ACGT]/) print header "\n" seq}' $OUTPUT_DIR/reads_amplicon_uniline.fasta > $OUTPUT_DIR/reads_amplicon_ACGT.fasta

		swarm \
			--differences 1 \
    			--fastidious \
    			--usearch-abundance \
    			--threads 8 \
    			--statistics-file $OUTPUT_DIR/SWARM.stats \
    			--internal-structure $OUTPUT_DIR/SWARM.struct \
    			--uclust-file $OUTPUT_DIR/SWARM.uc \
    			--seeds $OUTPUT_DIR/SWARM_representatives.fasta \
			$OUTPUT_DIR/reads_amplicon_ACGT.fasta \
			> $OUTPUT_DIR/SWARM.swarms

		rm $OUTPUT_DIR/reads_amplicon_uniline.fasta
	fi
fi

if [ "$CLUSTER" = "vsearch" ];then
	FINAL_FASTA="$OUTPUT_DIR/reads_OTU97_final.fasta"
elif [ "$CLUSTER" = "swarm" ];then
	FINAL_FASTA="$OUTPUT_DIR/SWARM_representatives.fasta"
fi

if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_OTU.uchime" ] &&  [ -s "$OUTPUT_DIR/reads_OTU_nonchimeras.fasta" ]; then
	echo "Output file already exists, skipping this step."
else
	echo "Running de novo chimera detection with vsearch on: $FINAL_FASTA"
	vsearch --uchime_denovo "$FINAL_FASTA" \
		--uchimeout $OUTPUT_DIR/reads_OTU.uchime \
		--nonchimeras $OUTPUT_DIR/reads_OTU_nonchimeras.fasta
fi

if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/taxonomy_OTU_usearch.txt" ] || [ -s "$OUTPUT_DIR/taxonomy_OTU_sintax.txt" ];then
	echo "Output file already exists, skipping this step."
else
	if [ "$METHOD" = "usearch" ]; then
		echo "Running usearch assignation with vsearch on: $OUTPUT_DIR/reads_OTU_nonchimeras.fasta"
		vsearch --usearch_global $OUTPUT_DIR/reads_OTU_nonchimeras.fasta \
    			--threads $nb_cores \
    			--dbmask none \
    			--qmask none \
    			--rowlen 0 \
    			--notrunclabels \
    			--userfields query+id+target \
    			--maxaccepts 10 \
    			--maxrejects 32 \
    			--db $path_to_db \
    			--id 0.7 \
    			--iddef 2 \
    			--userout $OUTPUT_DIR/taxonomy_OTU_usearch.txt

	elif [ "$METHOD" = "sintax" ]; then
		echo "Running sintax assignation with vsearch on: $OUTPUT_DIR/reads_OTU_nonchimeras.fasta"
		vsearch --sintax $OUTPUT_DIR/reads_OTU_nonchimeras.fasta   \
			--db $path_to_db  \
        		--sintax_cutoff $CUTOFF_PROB \
        		--tabbedout  $OUTPUT_DIR/taxonomy_OTU_sintax.txt    \
        		--strand plus    \
        		--threads $nb_cores
	fi
fi


if [ "$CLUSTER" = "vsearch" ];then
	STATS="$OUTPUT_DIR/stats_file_OTU97.txt"
	OTUS="$OUTPUT_DIR/reads_mapped_OTU97.txt"
	REPRESENTATIVES="$OUTPUT_DIR/reads_OTU97_final.fasta"
	UCHIME="$OUTPUT_DIR/reads_OTU.uchime"

elif [ "$CLUSTER" = "swarm" ];then
	REPRESENTATIVES="$OUTPUT_DIR/SWARM_representatives.fasta"
    	STATS="$OUTPUT_DIR/SWARM.stats"
    	OTUS="$OUTPUT_DIR/SWARM.swarms"
    	UCHIME="$OUTPUT_DIR/reads_OTU.uchime"
fi

if [ "$METHOD" = "sintax" ]; then
	ASSIGNMENTS="$OUTPUT_DIR/taxonomy_OTU_sintax.txt"
elif [ "$METHOD" = "usearch" ]; then
	cat $OUTPUT_DIR/taxonomy_OTU_usearch.txt | sort -k2 -n > $OUTPUT_DIR/taxonomy_OTU_usearch_sorted.txt
	ASSIGNMENTS="$OUTPUT_DIR/taxonomy_OTU_usearch_sorted.txt"
fi

OTU_TABLE="$OUTPUT_DIR/OTU_table_OTU.txt"

SCRIPT=$path_scripts"/OTU_contingency_table.py" 

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    "${METHOD}" \
    $INPUT_FILES/*.f* > "${OTU_TABLE}"
   
FILTERED="${OTU_TABLE/.txt/_filtered.txt}"
head -n 1 "${OTU_TABLE}" > "${FILTERED}"
cat "${OTU_TABLE}" | awk '$5 == "N" && $4 >= 100 && $2 > 1' >> "${FILTERED}"

rm $OUTPUT_DIR/taxonomy_OTU_usearch_sorted.txt
#python3 $DIR/$path_scripts/identity_distribution.py ../$OTU_TABLE

# cd $DIR
