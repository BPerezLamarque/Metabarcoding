#!/bin/bash



# This script performs OTU clustering using VSEARCH and generates a contingency table.
# It requires the following software: VSEARCH 2.29.3, Python 3.11.1.
# Depending on the options, it may also require swarm 3.1.3, SeqKit 2.9.0 and R 4.4.0. 
# It also requires the following parameters: path to the directory containing input fasta files, path to the database, path to the scripts, output directory, abundance threshold, and identity percentage.


NB_CORES=1
PATH_PYTHON=$(pwd)"/python_scripts/"
UNOISE="true"
CLUSTER="vsearch"
IDENTITY_PERCENTAGE=0.97
METHOD="vsearch"
CUTOFF_PROB=0.5
RESUME="true"
TAGJUMP="true"
MIN_READS=2

while getopts "hi:o:s:u:v:c:m:d:p:r:t:x:n:" option
do
        case $option in
                h)
                    echo "Usage: $0 -i <input_dir> -s <path_scripts> -d <path_db> -o <output_directory> -c <identity_percentage> -m <method> -p <proba_cutoff> -r <resume> -v <cluster_method> -t <tagjump>"
		    		echo "E.g: sbatch 02_OTU_clustering.sh -i Path_to/Demultiplexed_data/ -d Path_to_database/*.fasta -c 0.97 -o OTU_CLUSTERING -p 1 -u "false" -m "sintax" "
                    echo "  -i  Input directory containing FASTA files (can be compressed)"
                    echo "  -o  Output directory (default: 02_OTU_CLUSTERING)"
					echo "  -s  Path to the Python scripts"
		   		 	echo "  -u  Perform a denoising step with UNOISE before clustering (default = true)"
		    		echo "  -v  Clustering method: vsearch or swarm (default = "vsearch")"
		    		echo "  -c  Identity percentage for clustering with "vsearch" (default: 0.97)"
		    		echo "  -m  Method for taxonomic assignation: sintax or vsearch (default: "vsearch")"
		    		echo "  -d  Path to the database for taxonomic assignment"
		    		echo "  -p  SINTAX probability cutoff  (default: 0.5)"
		    		echo "  -r  Resume where the last output file was created (default: true)"
		    		echo "  -t  Perform a TagJump filtration (default = true)"
		    		echo "  -x  Minimum number of reads of an OTU to perform taxonomic assignment"
                    echo "  -h  Show this help message"		    		
                    echo "  -n  Number of CPU cores to use (default = 1)"
                    exit 0
                    ;;
                i)
		    		INPUT_DIR="$OPTARG"
                    ;;
                d)
                    path_to_db="$OPTARG"
                    ;;
                s)
                    PATH_PYTHON="$OPTARG"
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
		    	x)
		    		MIN_READS="$OPTARG"
		    		;;
		    	n)
		    		NB_CORES="$OPTARG"
		    		;;
        esac
done

if [ -z "$path_to_db" ]; then
    echo "Warning: No database provided. Skipping taxonomic assignment."
fi

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$(pwd)/02_OTU_CLUSTERING"
    mkdir -p "$OUTPUT_DIR"
fi

if [ "$CLUSTER" = "swarm" ];then
	UNOISE=false
fi


mkdir -p $OUTPUT_DIR


### DEREPLICATION

if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_amplicon_derep.fasta" ] && [ -s "$OUTPUT_DIR/reads_amplicon_derep.uc" ]; then
        echo "Output file already exists, skipping this step of dereplication."
else
	echo "Running Dereplication by vsearch on: $OUTPUT_DIR/reads_amplicon.fasta"
	cat $INPUT_DIR/*.f* > $OUTPUT_DIR/reads_amplicon.fasta
	vsearch \
		--derep_fulllength $OUTPUT_DIR/reads_amplicon.fasta \
    		--sizein \
    		--sizeout \
    		--relabel_sha1 \
    		--fasta_width 0 \
			--uc $OUTPUT_DIR/reads_amplicon_derep.uc \
    		--output $OUTPUT_DIR/reads_amplicon_derep.fasta
fi

### TAG-JUMP REMOVAL
# NEED dereplication.parquet and SeqTab.txt with sampleID seqID abundance

if [ "$TAGJUMP" = "true" ]; then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/Seq_tab_TagJumpFiltered.fasta" ]; then
                echo "Output file already exists, skipping this step of TagJump removal."
	else
		echo "Running the step of TagJump removal."
		for f in "$INPUT_DIR"/*.f*; do
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

		Rscript $PATH_PYTHON/tag_jump_removal_longtab.R \
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
		mv $OUTPUT_DIR/Seq_tab_TagJumpFiltered_derep.fasta $OUTPUT_DIR/reads_amplicon_derep.fasta
	fi
else
	 echo "Skipping the TagJump step."
fi


### DENOISING

if [ "$UNOISE" = "true" ]; then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/UNOISE.uc" ]; then
		echo "Output file already exists, skipping this step of denoising."
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
      			--threads $NB_CORES \
      			--fasta_width 0 \
      			--sizein --sizeout \
      			--centroids $OUTPUT_DIR/reads_amplicon_derep_clean.fasta \
      			--uc $OUTPUT_DIR/UNOISE.uc
      	
      	grep -v "^C" $OUTPUT_DIR/UNOISE.uc > $OUTPUT_DIR/UNOISE_hits.uc

	fi
else
	echo "Skipping the unoise step."
	cp $OUTPUT_DIR/reads_amplicon_derep.fasta $OUTPUT_DIR/reads_amplicon_derep_clean.fasta
fi


### SORT BY SIZE 

if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_amplicon_sorted.fasta" ]; then
	echo "Output file already exists, skipping this step of sorting by size."
else
	echo "Running sortbysize by vsearch on: $OUTPUT_DIR/reads_amplicon_derep_clean.fasta"
	vsearch -sortbysize $OUTPUT_DIR/reads_amplicon_derep_clean.fasta -output $OUTPUT_DIR/reads_amplicon_sorted.fasta -minsize 1
	
fi

### CLUSTERING

if [ "$CLUSTER" = "vsearch" ];then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_OTU_final.fasta" ] &&  [ -s "$OUTPUT_DIR/clusters_OTU.uc" ] && [ -s "$OUTPUT_DIR/reads_mapped_OTU.txt" ] &&  [ -s "$OUTPUT_DIR/stats_file_OTU.txt" ]; then
		echo "Output file already exists, skipping this step of clustering."
	else
		echo "Running vsearch clusterisation on: $OUTPUT_DIR/reads_amplicon_sorted.fasta"
		vsearch -cluster_size  $OUTPUT_DIR/reads_amplicon_sorted.fasta \
				--threads $NB_CORES \
    			--id $IDENTITY_PERCENTAGE \
    			--centroids $OUTPUT_DIR/reads_OTU.fasta \
    			--uc $OUTPUT_DIR/clusters_OTU.uc \
    			--sizein \
    			--sizeout
	
		vsearch --fasta_width 0 \
				--sortbysize $OUTPUT_DIR/reads_OTU.fasta \
				--output $OUTPUT_DIR/reads_OTU_final.fasta
				
				
		# Merge the two UC files (map all original reads to OTU clusters)
		if [ "$UNOISE" = "true" ];then
			#python3 $PATH_PYTHON/chain_uc.py $OUTPUT_DIR/UNOISE.uc $OUTPUT_DIR/clusters_OTU.uc $OUTPUT_DIR/clusters_OTU_full.uc
			python3 $PATH_PYTHON/chain_uc.py $OUTPUT_DIR/UNOISE_hits.uc $OUTPUT_DIR/clusters_OTU.uc $OUTPUT_DIR/clusters_OTU_full.uc
		else
			cat $OUTPUT_DIR/clusters_OTU.uc > $OUTPUT_DIR/clusters_OTU_full.uc
		fi


		python3 $PATH_PYTHON/map2qiime.py $OUTPUT_DIR/clusters_OTU_full.uc > $OUTPUT_DIR/reads_mapped_OTU.txt
		python3 $PATH_PYTHON/make_stats.py $OUTPUT_DIR/reads_OTU_final.fasta > $OUTPUT_DIR/stats_file_OTU.txt
	fi
elif [ "$CLUSTER" = "swarm" ];then
	if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/SWARM_representatives.fasta" ] &&  [ -s "$OUTPUT_DIR/SWARM.uc" ] && [ -s "$OUTPUT_DIR/SWARM.struct" ] &&  [ -s "$OUTPUT_DIR/SWARM.stats" ]; then
                echo "Output file already exists, skipping this step of clustering."
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
    			--threads $NB_CORES \
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
	FINAL_FASTA="$OUTPUT_DIR/reads_OTU_final.fasta"
elif [ "$CLUSTER" = "swarm" ];then
	FINAL_FASTA="$OUTPUT_DIR/SWARM_representatives.fasta"
fi


### DE NOVO CHIMERA FILTERING

if [ "$RESUME" = "true" ] && [ -s "$OUTPUT_DIR/reads_OTU.uchime" ] &&  [ -s "$OUTPUT_DIR/reads_OTU_nonchimeras.fasta" ]; then
	echo "Output file already exists, skipping this step of de novo chimera filtering."
else
	echo "Running de novo chimera detection with vsearch on: $FINAL_FASTA"
	vsearch --uchime_denovo "$FINAL_FASTA" \
		--uchimeout $OUTPUT_DIR/reads_OTU.uchime \
		--nonchimeras $OUTPUT_DIR/reads_OTU_nonchimeras.fasta
fi



# Only perform taxonomic assignation for OTU representing a certain number of reads
awk -v MIN_READS="$MIN_READS" '
  /^>/ {
    keep = 0
    if (match($0, /size=([0-9]+)/, a) && a[1] >= MIN_READS)
      keep = 1
  }
  keep
' "$OUTPUT_DIR/reads_OTU_nonchimeras.fasta" \
> "$OUTPUT_DIR/reads_OTU_nonchimeras_min_reads.fasta"


### TAXONOMIC ASSIGNATIONS

if [ "$RESUME" = "true" ] && \
   { [ -s "$OUTPUT_DIR/taxonomy_OTU_vsearch_sorted.txt" ] || \
     [ -s "$OUTPUT_DIR/taxonomy_OTU_sintax.txt" ]; }; then
     	echo "Output file already exists, skipping this step of taxonomic assignations."
else
	if [ "$METHOD" = "vsearch" ]; then
		echo "Running vsearch assignation with vsearch on: $OUTPUT_DIR/reads_OTU_nonchimeras_min_reads.fasta"
		vsearch --usearch_global $OUTPUT_DIR/reads_OTU_nonchimeras_min_reads.fasta \
    			--threads $NB_CORES \
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
    			--userout $OUTPUT_DIR/taxonomy_OTU_vsearch.txt

	elif [ "$METHOD" = "sintax" ]; then
		echo "Running sintax assignation with vsearch on: $OUTPUT_DIR/reads_OTU_nonchimeras_min_reads.fasta"
		vsearch --sintax $OUTPUT_DIR/reads_OTU_nonchimeras_min_reads.fasta   \
			--db $path_to_db  \
        		--sintax_cutoff $CUTOFF_PROB \
        		--tabbedout  $OUTPUT_DIR/taxonomy_OTU_sintax.txt    \
        		--strand plus    \
        		--threads $NB_CORES
	fi
fi



### PREPARE THE OTU TABLE 

if [ "$CLUSTER" = "vsearch" ];then
	STATS="$OUTPUT_DIR/stats_file_OTU.txt"
	OTUS="$OUTPUT_DIR/reads_mapped_OTU.txt"
	REPRESENTATIVES="$OUTPUT_DIR/reads_OTU_final.fasta"
	UCHIME="$OUTPUT_DIR/reads_OTU.uchime"

elif [ "$CLUSTER" = "swarm" ];then
	REPRESENTATIVES="$OUTPUT_DIR/SWARM_representatives.fasta"
    STATS="$OUTPUT_DIR/SWARM.stats"
    OTUS="$OUTPUT_DIR/SWARM.swarms"
    UCHIME="$OUTPUT_DIR/reads_OTU.uchime"
fi

if [ "$METHOD" = "sintax" ]; then
	ASSIGNMENTS="$OUTPUT_DIR/taxonomy_OTU_sintax.txt"
elif [ "$METHOD" = "vsearch" ]; then
	# Keep one line per OTU
	awk '{
		  otu=$1; id=$2
		  if (!(otu in best) || id>best[otu]) { best[otu]=id; line[otu]=$0 }
		  if (!(otu in seen_order)) { seen_order[otu]=NR; order[NR]=otu }
	} END {
		  for(i=1;i<=NR;i++){ otu=order[i]; if(otu in line){ print line[otu]; delete line[otu] } }
	}' "$OUTPUT_DIR/taxonomy_OTU_vsearch.txt" > "$OUTPUT_DIR/taxonomy_OTU_vsearch_sorted.txt"
	ASSIGNMENTS="$OUTPUT_DIR/taxonomy_OTU_vsearch_sorted.txt"
fi


### MAKE THE OTU TABLE 

OTU_TABLE="$OUTPUT_DIR/OTU_table.txt"

SCRIPT=$PATH_PYTHON"/OTU_contingency_table.py" 

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    "${METHOD}" \
    $INPUT_DIR/*.f* > "${OTU_TABLE}"
   
FILTERED="${OTU_TABLE/.txt/_filtered.txt}"
head -n 1 "${OTU_TABLE}" > "${FILTERED}"
cat "${OTU_TABLE}" | awk '$5 == "N" && $4 >= 50' >> "${FILTERED}"

#In order to prepare the LULU step
cut -f3,9- "${FILTERED}" > $OUTPUT_DIR/pre_LULU_match.txt


# remove the intermediary files: 

remove=true

if [ "$remove" = true ]; then

    rm -f "$OUTPUT_DIR"/reads_amplicon*
    rm -f "$OUTPUT_DIR"/Seq*
    rm -f "$OUTPUT_DIR/TagJump_scores.rds"
    rm -f "$OUTPUT_DIR/reads_OTU.uchime"
    rm -f "$OUTPUT_DIR/reads_OTU_nonchimeras_min_reads.fasta"
    rm -f "$FINAL_FASTA"

    if [ "$CLUSTER" = "vsearch" ]; then
        rm -f "$OUTPUT_DIR/UNOISE.uc"
        rm -f "$OUTPUT_DIR/UNOISE_hits.uc"
    	rm -f "$OUTPUT_DIR/reads_OTU.fasta"
        rm -f "$OUTPUT_DIR/clusters_OTU.uc"
        rm -f "$OUTPUT_DIR/clusters_OTU_full.uc.uc"
    fi

    if [ "$CLUSTER" = "swarm" ]; then
        rm -f "$OUTPUT_DIR/SWARM.uc"
        rm -f "$OUTPUT_DIR/SWARM.struct"
        rm -f "$OUTPUT_DIR/SWARM.stats"
        rm -f "$OUTPUT_DIR/SWARM.swarms"
    fi

    rm -f "$OUTPUT_DIR/taxonomy_OTU_vsearch.txt"
fi

