# Metabarcoding


This GitHub repository provides pipelines for conducting bioinformatic analyses of fungal metabarcoding using both Illumina short reads and PacBio long reads. The Illumina workflow is derived from the methods developed by Frédéric Mahé (https://github.com/torognes/vsearch), whereas the PacBio workflow (NextITS) is based on the work of Vladimir Mikryukov (https://github.com/vmikk/NextITS; documentation: https://next-its.github.io/).
The development of these pipelines was supported by funding from the MITI DEFIS project (PI: Mélanie Roy).

**Contributors:** Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com), Valentin Etienne, and Florian Tilliet. 





# Step 0: Preparing the database for taxonomic assignment (Illumina and PacBio)

This step prepares a reference database for taxonomic assignment by trimming primer sequences from the input FASTA file using **Cutadapt**.

**Requirements:**
* **Software:** Cutadapt ≥ 5.0
* **Inputs:**
  * Raw database FASTA file (compressed or uncompressed)
  * Forward primer sequence
  * Reverse primer sequence

### Running the script:

```bash
INPUT="SINTAX_EUK_ITS_v1.9.4.fasta"
PATH_DIR_DB="Path_to/database/"
NAME_DB="General_EUK_ITS_v1.9.4_Ted"

bash cut_database.sh \
    -i "$INPUT" \
    -f GGAAGTAAAAGTCGTAACAAGG \
    -r CAAGAGATCCGTTGTTGAAAGTK \
    -o "$PATH_DIR_DB" \
    -n "$NAME_DB" \
    -x false \
    -y false \
    -l 100 \
    -t "sintax"
```




# Analyzing an Illumina dataset

## Step 1: Demultiplexing and quality control

This step merges paired-end reads, demultiplexes them, and dereplicates sequences at the sample level. It also performs quality filtering and evaluates the quality of the merged reads using **FastQC**.

**Required software:**
* Cutadapt ≥ 5.0
* VSEARCH ≥ 2.29.3
* FastQC ≥ 0.12.1

**Required inputs:**
* Path to the directory containing the raw FASTQ files
* A mapping file for demultiplexing (optional)

If demultiplexing is needed, the mapping file (`mapping_file.txt`) must contain the following columns:
* `SampleID`
* `barcodeFw`
* `primerFw`
* `primerRev`
* `barcodeRev`

**Format example:**

```
Sample   barcodeRev   barcodeFw   primerRev                   primerFw                   combineFw                                  combineRev
B01      ACACACAC     ACACACAC    CAAGAGATCCGTTGTTGAAAGTK     GGAAGTAAAAGTCGTAACAAGG     ACACACACCAAGAGATCCGTTGTTGAAAGTK            GGAAGTAAAAGTCGTAACAAGGACACACAC
M5_10    ACAGCACA     ACACACAC    CAAGAGATCCGTTGTTGAAAGTK     GGAAGTAAAAGTCGTAACAAGG     ACAGCACACAAGAGATCCGTTGTTGAAAGTK            GGAAGTAAAAGTCGTAACAAGGACACACAC
M7_4     GTGTACAT     ACACACAC    CAAGAGATCCGTTGTTGAAAGTK     GGAAGTAAAAGTCGTAACAAGG     GTGTACATCAAGAGATCCGTTGTTGAAAGTK            GGAAGTAAAAGTCGTAACAAGGACACACAC
```

If no mapping file is provided, the demultiplexing step will be skipped.

### Running the script:

```bash
INPUT_FILES="path/to/directory/*.fastq"  # Corresponds to the Illumina R1 and R2 FASTQ files

bash 01_demultiplexing_QC.sh \
    -i "$INPUT_FILES" \
    -o ./01_OTU_PRECLUSTERING \
    -m mapping_file.txt
```





## Step 2: OTU clustering and taxonomic assignment

This step performs OTU clustering using **VSEARCH** and generates a OTU table. Depending on the selected options, it can also run **UNOISE denoising**, **swarm clustering**, and **taxonomic assignment** using **SINTAX** or **VSEARCH**.

**Required software:**
* VSEARCH ≥ 2.29.3
* Python ≥ 3.11.1

**Additional software depending on options:**
* swarm ≥ 3.1.3
* SeqKit ≥ 2.9.0
* R ≥ 4.4.0

**Required inputs and parameters:**

```
-i  Input directory containing FASTA files (compressed or uncompressed)
-o  Output directory (default: 02_OTU_CLUSTERING)
-s  Path to the Python scripts
-u  Perform a denoising step with UNOISE before clustering (default: true)
-v  Clustering method: "vsearch" or "swarm" (default: "vsearch")
-c  Identity percentage for clustering (default: 0.97)
-m  Method for taxonomic assignment: "sintax" or "vsearch" (default: "vsearch")
-d  Path to the database for taxonomic assignment
-p  SINTAX probability cutoff (default: 0.5)
-n  Number of CPU cores to use (default: 1)
```

### Run the script:

Example of script using VSEARCH OTU clustering at 97% identity, with UNOISE denoising, and SINTAX taxonomic assignment (probability cutoff = 0.5):

```bash
OUT_DIR="02_OTU_CLUSTERING_VSEARCH_sintax"

bash 02_OTU_clustering \
    -i 01_OTU_PRECLUSTERING/Demultiplexed_data/ \
    -s Path_to/python_scripts/ \
    -u "true" \
    -v "vsearch" \
    -c 0.97 \
    -m "sintax" \
    -d "$PATH_DIR_DB/$NAME_DB" \
    -p 0.5 \
    -n 1 \
    -o "$OUT_DIR"
```





## Step 3 — Running LULU for OTU table curation

This step refines the OTU table using **LULU**, which relies on sequence similarity and co-occurrence patterns to curate inflated OTU sets.

**Required software:**
* VSEARCH ≥ 2.29
* gcc ≥ 10
* **mumu** ([https://github.com/frederic-mahe/mumu.git](https://github.com/frederic-mahe/mumu.git))

**Installing mumu:**
```bash
git clone https://github.com/frederic-mahe/mumu.git
cd mumu
module load compilers/gcc/12.2.0
make
```

You can verify the installation with:
```bash
./mumu --help
```

### Running the script:

```bash
bash 03_LULU_curation.sh \
    -f "$OUT_DIR/reads_OTU_nonchimeras.fasta" \
    -t "$OUT_DIR/pre_LULU_match.txt" \
    -o SWARM_USEARCH_LULU_OTU_CLUSTERING \
    -c "$OUT_DIR/OTU_table_OTU_filtered.txt"
```




# Analyzing a PacBio dataset: 

The scripts 1 to 3 can be run successively for 




    
    