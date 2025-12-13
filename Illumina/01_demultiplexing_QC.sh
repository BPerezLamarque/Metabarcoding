#!/bin/bash

## This script merges paired-end reads, demultiplexes them, and dereplicates the sequences at the sample level. It also performs quality filtering and checks the quality of the merged reads using FastQC.
## The following software is required: Cutadapt 5.0, VSEARCH 2.29.3, FastQC 0.12.1.
## It requires the following parameters: path to the directory containing input fastq files, and mapping file.
## For demultiplexing, it requires a mapping file with the following columns: SampleID, barcodeFw, primerFw, primerRev, barcodeRev.



while getopts "o:i:d:m:f:r:" option
do
        case $option in
                h)
                    echo "Usage: $0 -i <input_files> -o <output_directory> -m <mapping_file>"
                    echo "  -i  Input directory containing fastq files (can be compressed; can be a part of the files names, e.g. '../raw_reads/ITS*.fastq')"
                    echo "  -o  Output directory (default: ./01_OTU_PRECLUSTERING)"
                    echo "  -m  Mapping file (tab-separated) with the following columns: SampleID, barcodeFw, primerFw, primerRev, barcodeRev"
                    echo "  -h  Show this help message"
                    echo "  -l  minimum length of the retained sequences"
                    echo "  -n  number of cores"
                    exit 0
                    ;;
                o)
                    OUTPUT_DIR="$OPTARG"
                    ;;
                i)
		    		INPUT_FILES="$OPTARG"
                    ;;
                m)
                    MAPPING="$OPTARG"
		    		;;
		    	l)
                    MIN_LENGTH="$OPTARG"
                    ;;
                n)
                    NB_CORES="$OPTARG"
                    ;;
        esac
done


if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$(pwd)/01_OTU_PRECLUSTERING"
fi

if [ -z "$MIN_LENGTH" ]; then
    MIN_LENGTH=80
fi

if [ -z "$NB_CORES" ]; then
    NB_CORES=1
fi

# Prepare the data
echo "Prepare the data"

mkdir -p $OUTPUT_DIR

# List all samples, remove the extension, and store the names
ls $INPUT_FILES | sed 's#.*/##' | sed 's/.f.*//' | sed 's/_R[12].*//' | sort -u > $OUTPUT_DIR/list_sample.txt

# Check the quality encoding (33 or 64?)
OUTPUT="$OUTPUT_DIR/Quality.encoding.log"
FIRST_SAMPLE=$(head -n 1 $OUTPUT_DIR/list_sample.txt)
INPUT_SAMPLE=$(ls $INPUT_FILES | grep -E "${FIRST_SAMPLE}_R1\.(fastq\.gz|fastq)$")

if file "$INPUT_SAMPLE" | grep -q "compressed"; then
    zcat "$INPUT_SAMPLE" > $OUTPUT_DIR/temp_decompressed.fastq
else
    cat "$INPUT_SAMPLE" > $OUTPUT_DIR/temp_decompressed.fastq
fi
vsearch \
    --threads $NB_CORES \
    --fastq_chars $OUTPUT_DIR/temp_decompressed.fastq 2> ${OUTPUT}
QUALITY_ENCODING=$(grep "Guess: Original" $OUTPUT | awk -F'[+]' '{print $2}' | awk -F')' '{print $1}')


# Merge the paired-end reads
echo "Merge the paired-end reads"

mkdir -p $OUTPUT_DIR/merged_reads/

for sample in $(cat $OUTPUT_DIR/list_sample.txt); do
    
    INPUT_R1=$(ls $INPUT_FILES | grep -E "${sample}_R1\.(fastq\.gz|fastq)$")
    INPUT_R2=$(ls $INPUT_FILES | grep -E "${sample}_R2\.(fastq\.gz|fastq)$")

    OUTPUT="$OUTPUT_DIR/merged_reads/"$sample".fastq"

    vsearch \
        --fastq_mergepairs ${INPUT_R1} \
        --reverse ${INPUT_R2} \
        --fastq_allowmergestagger \
        --fastq_ascii $QUALITY_ENCODING \
        --fastqout $OUTPUT \
        --quiet 2>> ${OUTPUT/.fastq/.log}
    
done


# Checking the quality with FASTQC
echo "Checking the quality with FASTQC"

# FastQC is available on https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip)

mkdir -p $OUTPUT_DIR/FastQC/

for sample in $(cat $OUTPUT_DIR/list_sample.txt); do
    
    INPUT="$OUTPUT_DIR/merged_reads/"$sample".fastq"
    
    fastqc ${INPUT} -t 8 -o $OUTPUT_DIR/FastQC/
    
done

# Quality filtering
echo "Quality filtering"

for sample in $(cat $OUTPUT_DIR/list_sample.txt); do

    echo $sample

    INPUT="$OUTPUT_DIR/merged_reads/"$sample".fastq"   # fastq
    OUTPUT="$OUTPUT_DIR/merged_reads/"$sample".fasta"  # fasta

    vsearch \
    --fastq_filter $INPUT \
    --fastq_maxns 0 \
    --fastq_maxee 2 \
    --fastaout "${OUTPUT}"
    
    echo "number of reads in $sample before quality filtering:" > "$OUTPUT_DIR/merged_reads/"$sample"_QF.log" 
    cat $INPUT | wc -l | awk '{print int($1/4)}' >> "$OUTPUT_DIR/merged_reads/"$sample"_QF.log" 
    echo "number of reads in $sample after quality filtering:" >> "$OUTPUT_DIR/merged_reads/"$sample"_QF.log" 
    grep ">" "${OUTPUT}" | wc -l >> "$OUTPUT_DIR/merged_reads/"$sample"_QF.log" 

done

# Check if the mapping file is specified, if not, dereplicate the sequences at the sample level and end here
if [ -z "$MAPPING" ]; then
    mkdir -p $OUTPUT_DIR/merged_derep_reads/
    for sample in $(cat $OUTPUT_DIR/list_sample.txt); do

    echo $sample

    INPUT="$OUTPUT_DIR/merged_reads/"$sample".fasta"  # fasta
    OUTPUT="$OUTPUT_DIR/merged_derep_reads/"$sample"_dereplicated.fasta"  # fasta

     # Dereplicate at the sample level
    vsearch --quiet \
        --threads $NB_CORES \
        --derep_fulllength "temp_"$sample".fasta" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --output "${OUTPUT}"

    done
    exit 0
fi

# Demultiplexing
echo "Demultiplexing"

# Get the column index of the barcode and primer columns

barcodeFwColumnIdx="$(( $(head -n 1 $MAPPING | tr '\t' '\n' | grep -nx 'barcodeFw' | cut -d: -f1) -1 ))" # getting the index of column with forward barcode etc. "-1" is applied as array below counts from 0
primerFwColumnIdx="$(( $(head -n 1 $MAPPING | tr '\t' '\n' | grep -nx 'primerFw' | cut -d: -f1) -1 ))"
barcodeRevColumnsIdx="$(( $(head -n 1 $MAPPING | tr '\t' '\n' | grep -nx 'barcodeRev' | cut -d: -f1) -1 ))"
primerRevColumnsIdx="$(( $(head -n 1 $MAPPING | tr '\t' '\n' | grep -nx 'primerRev' | cut -d: -f1) -1 ))"
 
# Create the output directory 
mkdir -p "$OUTPUT_DIR/Demultiplexed_data/"
mkdir -p "$OUTPUT_DIR/tmp_demux"


INPUT=$(find "$OUTPUT_DIR/merged_reads/" -type f -name "*.fasta" ! -name "*_RC.fasta" | head -n 1)
INPUT_REVCOMP="${INPUT/.fasta/_RC.fasta}"


# Reverse complement fastq file
vsearch --quiet \
    --threads $NB_CORES \
    --fastx_revcomp "${INPUT}" \
    --fastaout "${INPUT_REVCOMP}"

while read -r line; do
    if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
        IFS=$'\t' read -r -a array <<< "$line"

        # Get sequences
        FwBarcode="${array[${barcodeFwColumnIdx}]}"
        FwPrimer="${array[${primerFwColumnIdx}]}"
        RevBarcode="${array[${barcodeRevColumnsIdx}]}"
        RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
        RevPrimer="${array[${primerRevColumnsIdx}]}"
        RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
        SAMPLE_NAME="${array[0]}"
        
        # Output file names
        LOG="$OUTPUT_DIR/Demultiplexed_data/${SAMPLE_NAME}.log"
        FINAL_FASTA="$OUTPUT_DIR/Demultiplexed_data/${SAMPLE_NAME}.fas"
        
        # Some information
        echo "${SAMPLE_NAME} is being processed.."
        echo "Barcode Fw: ${FwBarcode}"
        echo "Primer Fw: ${FwPrimer}"
        echo "Primer Rev (RC): ${RevPrimerRC}"
        echo "Barcode Rev (RC): ${RevBarcodeRC}"

        if [ -f "${FINAL_FASTA}" ]; then
            echo "${SAMPLE_NAME} has already been processed. Skipping..."
            continue
        fi
        
        function trim_without_ambiguity {

            SEQTOT="${FwBarcode}${FwPrimer}${RevPrimerRC}${RevBarcodeRC}"
            MIN_MATCHED=${#SEQTOT}
            ERROR_RATE=0

            cat "${INPUT}" "${INPUT_REVCOMP}" | cutadapt --cores 0 -g "${FwBarcode}${FwPrimer}...${RevPrimerRC}${RevBarcodeRC}" --discard-untrimmed --minimum-length "${MIN_LENGTH}" -O ${MIN_MATCHED} -e "${ERROR_RATE}" - 2> "${LOG}" > "$OUTPUT_DIR/tmp_demux/temp_${SAMPLE_NAME}.fasta"
        }

        trim_without_ambiguity

        # Dereplicate at the study level
        vsearch --quiet \
            --threads $NB_CORES \
            --derep_fulllength "$OUTPUT_DIR/tmp_demux/temp_${SAMPLE_NAME}.fasta" \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --relabel_sha1 \
            --output "${FINAL_FASTA}" 2>> "${LOG}"

    fi
done < "${MAPPING}"

# Note: the option sha1 (encoding system) is giving the same names to the identical amplicons across samples

# Remove the files that are no longer useful
rm -f "${INPUT}" "${INPUT_REVCOMP}" "temp.fastq" "temp.fasta" "${OUTPUT_DIR}/temp_decompressed.fastq" 
rm -rf "${OUTPUT_DIR}/tmp_demux/"
rm -f "${OUTPUT_DIR}/merged_reads/*fastq"
rm -f "${OUTPUT_DIR}/merged_reads/*fasta"
