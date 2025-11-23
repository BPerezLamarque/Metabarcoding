#!/bin/bash

# This script prepares a database for taxonomic assignation by trimming primers from the input sequences using Cutadapt.
# It requires the following software: Cutadapt 5.0.
# It also requires the following parameters: raw database input fasta file (can be compressed), forward primer, reverse primer.


while getopts "hi:o:f:r:n:x:y:l:" option
do
        case $option in
                h)
                    echo "Usage: $0 -i <input_file> -f <forward_primer> -r <reverse_primer> -o <output_directory> -n -x -y -l -t <type>"
                    echo "  -i  Input raw database file (FASTA format, can be compressed)"
                    echo "  -f  Forward primer sequence"
                    echo "  -r  Reverse primer sequence"
                    echo "  -o  Output directory (default: current directory)"
                    echo "  -n  Name of the database"
                    echo "  -x  indicates whether finding (and extracting) the forward primer is mandatory for keeping the sequence (true or false)"
                    echo "  -y  indicates whether finding (and extracting) the reverse primer is mandatory for keeping the sequence (true or false)"
                    echo "  -l  minimum length of the retained sequences"
		    echo "  -t  Database type: 'usearch' or 'sintax' (default: 'usearch')"
		    exit 0
                    ;;
                i)
                    INPUT="$OPTARG"
                    ;;
                f)
                    PRIMER_F="$OPTARG"
                    ;;
                r)
                    PRIMER_R="$OPTARG"
                    ;;
                o)
                    DIR="$OPTARG"
                    ;;
                n)
                    NAME="$OPTARG"
                    ;;
                x)
                    OBLIGATE_F="$OPTARG"
                    ;;
                y)
                    OBLIGATE_R="$OPTARG"
                    ;;
                l)
                    MIN_LENGTH="$OPTARG"
                    ;;
		t)
		    DB_TYPE="$OPTARG" 
		    ;;
        esac
done


if [ -z "$INPUT" ]; then
    echo "Error: An input file must be specified."
    exit 1
fi

if (file $INPUT | grep -q "compressed data");
then
    command="zcat"
else
    command="cat"
fi

if [ -z "$NAME" ]; then
    NAME=$INPUT
fi

if [ -z "$OBLIGATE_F" ]; then
    OBLIGATE_F=false
fi

if [ -z "$OBLIGATE_R" ]; then
    OBLIGATE_R=false
fi

if [ -z "$MIN_LENGTH" ]; then
    MIN_LENGTH=100
fi

if [ -z "$DIR" ]; then
    DIR="$(pwd)"
fi

if [ -t "$DB_TYPE" ]; then
    DB_TYPE="usearch"
fi

mkdir -p "$DIR"
cd $DIR


# Define variables and output files primers
OUTPUT=$DIR/${NAME%.*}_trimmed_temp.fasta
LOG=$DIR/${NAME%.*}_trimmed_temp.log
OUTPUT_R=$DIR/${NAME%.*}_trimmed.fasta
LOG_R=$DIR/${NAME%.*}_trimmed.log


# Reverse complement of the reverse primer
ANTI_PRIMER_R=$( echo "${PRIMER_R}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )


# Paramerers of cutadapt
MIN_F=$(( ${#PRIMER_F} * 2 / 3 )) # should match at least 66% of the forward primer
MIN_R=$(( ${#ANTI_PRIMER_R} * 2 / 3 )) # should match at least 66% of the reverse primer
ERROR_RATE=0.10


# Trim forward
if [ "$OBLIGATE_F" = true ]; then
    CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH} --no-indels -e ${ERROR_RATE}"
else
    CUTADAPT="cutadapt --minimum-length ${MIN_LENGTH} --no-indels -e ${ERROR_RATE}"
fi

if  [ "$DB_TYPE" = "usearch" ]; then
	${command} "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
		${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
		sed '/^>/ s/;/|/g ; /^>/ s/ /_/g' > "${OUTPUT}"
else
	${command} "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
        ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
		sed '/^>/ s/|/;/g ; /^>/ s/ /_/g' > "${OUTPUT}"


# Trim reverse
if [ "$OBLIGATE_R" = true ]; then
    CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH} --no-indels -e ${ERROR_RATE}"
else
    CUTADAPT="cutadapt --minimum-length ${MIN_LENGTH} --no-indels -e ${ERROR_RATE}"
fi

if  [ "$DB_TYPE" = "usearch" ]; then
	${command} "${OUTPUT}" | sed '/^>/ ! s/U/T/g' | \
		${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2> "${LOG_R}" | \
		sed '/^>/ s/;/|/g ; /^>/ s/ /_/g' > "${OUTPUT_R}"
else
	${command} "${OUTPUT}" | sed '/^>/ ! s/U/T/g' | \
     	${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2> "${LOG_R}" | \
		sed '/^>/ s/|/;/g ; /^>/ s/ /_/g' > "${OUTPUT_R}"
rm $OUTPUT
