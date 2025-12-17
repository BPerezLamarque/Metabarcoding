# Metabarcoding


This GitHub repository provides pipelines for conducting bioinformatic analyses of fungal metabarcoding using both Illumina short reads and PacBio long reads. The Illumina workflow is derived from the methods developed by Frédéric Mahé (https://github.com/torognes/vsearch), whereas the PacBio workflow (NextITS) is based on the work of Vladimir Mikryukov (https://github.com/vmikk/NextITS; documentation: https://next-its.github.io/).


The development of these pipelines was supported by funding from the MITI DEFIS project (PI: Mélanie Roy).


**Contributors:** Valentin Etienne, Florian Tilliet, and Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com). 



# Table of Contents

* [Step 0: Preparing the database](#step-0--preparing-the-database-for-taxonomic-assignment)
  
* [Analyzing an Illumina dataset](#analyzing-an-illumina-dataset)
  * [Step 1: Demultiplexing and quality control](#step-1:-demultiplexing-and-quality-control)
  * [Step 2: OTU clustering and taxonomic assignment](#step-2:-otu-clustering-and-taxonomic-assignment)
  * [Step 3: Running LULU for OTU table curation](#step-3:-running-lulu-for-otu-table-curation)

* [Analyzing a PacBio dataset](#analyzing-a-pacbio-dataset)
  * [Step 1: Quality control and artefact removal](#step-1:-quality-control-and-artefact-removal)
  * [Step 2: OTU clustering](#step-2:-otu-clustering-1)
  * [Step 3: Taxonomic assignation](#step-3:-taxonomic-assignation)




# Step 0: Preparing the database for taxonomic assignment (Illumina and PacBio):

This step (https://github.com/BPerezLamarque/Metabarcoding/blob/main/Illumina/00_prepare_database.sh) prepares a reference database for taxonomic assignment by trimming primer sequences from the input FASTA file using **Cutadapt**.



**Requirements:**
* **Software:** Cutadapt ≥ 5.0
* **Inputs:**
  * **-i**: Raw database FASTA file (compressed or uncompressed)
  * **-f**: Forward primer sequence
  * **-r**: Reverse primer sequence
  * **-x**: Whether detecting and trimming the forward primer is mandatory to keep the sequence (`true` or `false`)
  * **-y**: Whether detecting and trimming the reverse primer is mandatory to keep the sequence (`true` or `false`)
  * **-l**: Minimum sequence length required to retain the entry in the database

### Running the script:

```bash
INPUT="General_EUK_ITS_v2.0.fasta"  # Downloaded from EUKARYOME (https://eukaryome.org/generalfasta/)
NAME_DB="EUK_ITS1_fung02"
PATH_DIR_DB="Path_to/database/"


bash 00_prepare_database.sh \
    -i "$INPUT" \
    -f GGAAGTAAAAGTCGTAACAAGG \
    -r CAAGAGATCCGTTGTTGAAAGTK \
    -o "$PATH_DIR_DB" \
    -n "$NAME_DB" \
    -x false \
    -y false \
    -l 100
```

---

If you want, I can automatically harmonize all requirement sections across Steps 0, 1, 2, and PacBio to ensure identical structure.





# Analyzing an Illumina dataset:

## Step 1: Demultiplexing and quality control:

This step (https://github.com/BPerezLamarque/Metabarcoding/blob/main/Illumina/01_demultiplexing_QC.sh) merges paired-end reads, demultiplexes them, and dereplicates sequences at the sample level. It also performs quality filtering and evaluates the quality of the merged reads using **FastQC**.

**Required software:**
* Cutadapt ≥ 5.0
* VSEARCH ≥ 2.29.3
* FastQC ≥ 0.12.1

**Required inputs:**
* Path to the directory containing the raw FASTQ files
* A mapping file for demultiplexing (optional)
* A minimum sequence length required to keep the sequences during demultiplexing (**-l**)

If demultiplexing is needed, the mapping file (`mapping_file.txt`) must contain the following columns:
* `SampleID`
* `barcodeFw`
* `primerFw`
* `primerRev`
* `barcodeRev`

**Format example:**

```
Sample   barcodeFw    barcodeRev  primerFw                   primerRev                  combineFw                                  combineRev
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
    -m mapping_file.txt \
    -l 80
```





## Step 2: OTU clustering and taxonomic assignment:

This step (https://github.com/BPerezLamarque/Metabarcoding/blob/main/Illumina/02_OTU_clustering.sh) performs OTU clustering using **VSEARCH** and generates a OTU table. Depending on the selected options, it can also run **UNOISE denoising**, **swarm clustering**, and **taxonomic assignment** using **VSEARCH** or **SINTAX**  (NB: SINTAX requires a specific formating of the taxonomic ranks, following the taxonomy in UNITE).

**Required software:**
* VSEARCH ≥ 2.29.3
* Python ≥ 3.11.1

**Additional software depending on options:**
* swarm ≥ 3.1.3


**Required inputs and parameters:**

```
-i  Input directory containing FASTA files (compressed or uncompressed)
-o  Output directory (default: 02_OTU_CLUSTERING)
-s  Path to the Python scripts
-u  Perform a denoising step with UNOISE before clustering (default: true)
-v  Clustering method: "vsearch" or "swarm" (default: "vsearch")
-c  Identity percentage for clustering with "vsearch" (default: 0.97)
-m  Method for taxonomic assignment: "vsearch" or "sintax" (default: "vsearch")
-d  Path to the database for taxonomic assignment
-p  SINTAX probability cutoff (default: 0.5)
-x  Minimum number of reads of an OTU to perform taxonomic assignment
-n  Number of CPU cores to use (default: 1)
```


### Run the script:

Example of script using VSEARCH OTU clustering at 97% identity, with UNOISE denoising, and SINTAX taxonomic assignment (probability cutoff = 0.5):

```bash
OUT_DIR="02_OTU_CLUSTERING_OTU97_VSEARCH"

bash 02_OTU_clustering.sh \
    -i 01_OTU_PRECLUSTERING/Demultiplexed_data/ \
    -s Path_to/python_scripts/ \
    -u "true" \
    -v "vsearch" \
    -c 0.97 \
    -m "vsearch" \
    -d "$PATH_DIR_DB/$NAME_DB" \
    -x 2 \
    -n 1 \
    -o "$OUT_DIR"
```





## Step 3: Running LULU for OTU table curation:

This step (https://github.com/BPerezLamarque/Metabarcoding/blob/main/Illumina/03_LULU_curation.sh) refines the OTU table using **LULU**, which relies on sequence similarity and co-occurrence patterns to curate inflated OTU sets.

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
    -c "$OUT_DIR/OTU_table_OTU_filtered.txt" \
	-s 0.97
```






# Analyzing a PacBio dataset:

The PacBio workflow (NextITS) is based on the work of **Vladimir Mikryukov** ([https://github.com/vmikk/NextITS](https://github.com/vmikk/NextITS)).
Additional documentation for Steps 1 and 2 is available at: [https://next-its.github.io/](https://next-its.github.io/)

**Citation:**
Mikryukov V., Anslan S., Tedersoo L. *NextITS: a pipeline for metabarcoding fungi and other eukaryotes with full-length ITS sequenced with PacBio.*
[https://github.com/vmikk/NextITS](https://github.com/vmikk/NextITS) — DOI: 10.5281/zenodo.15074881



## Step 1: Quality control and artefact removal:

This step runs the first module of the NextITS pipeline, performing quality filtering, primer trimming, homopolymer correction, and artefact removal.

**Required software/modules:**

* Nextflow ≥ 25.04.0
* Singularity ≥ 3.9.9
* Java ≥ 17.0.6

### Getting and preparing NextITS

```bash
# Download the latest version of NextITS
nextflow pull vmikk/NextITS

# Download the container image containing all dependencies
singularity pull nextits.sif docker://vmikk/nextits:1.1.0

# Display available options and help
nextflow run vmikk/NextITS --help
```

### Running Step 1:

Adjust primer sequences and the `its_region` parameter if different from the defaults using the "full" ITS region.
Using a chimera reference database is optional (denovo detection only is also possible).

```bash
nextflow run vmikk/NextITS -r main \
  -resume \
  --demultiplexed true \
  --step "Step1" \
  --input "$(pwd)/path/to/fastq_files" \
  --primer_forward GTACACACCGCCCGTCG \
  --primer_reverse CGCCTSCSCTTANTDATATGC \
  --its_region "full" \
  --hp "true" \
  --chimera_methods "ref,denovo" \
  --chimera_db "$(pwd)/path/to/chimera_database" \
  --tj "true" \
  --outdir "Step1_Results/test" \
  -with-singularity "$(pwd)/nextits.sif"
```

---

## Step 2: OTU clustering:

This step performs denoising (UNOISE), OTU clustering at 97% identity using VSEARCH, and OTU curation using **LULU**.

```bash
nextflow run vmikk/NextITS -r main \
  --step "Step2" \
  --data_path "$(pwd)/Step1_Results/" \
  --preclustering "unoise" \
  --clustering "vsearch" \
  --otu_id 0.97 \
  --otu_iddef 2 \
  --lulu "true" \
  --outdir "Step2_vsearch" \
  -with-singularity "$(pwd)/../nextits.sif"
```

---

## Step 3: Taxonomic assignation:

This step performs taxonomic assignment using **VSEARCH**, with either the **SINTAX** or **usearch_global** algorithm.

**Required software:**
* VSEARCH ≥ 2.29


### Assignation with SINTAX:

```bash
vsearch --sintax Path/to/fasta_file_from_step_2 \
    --db "$PATH_DIR_DB/$NAME_DB" \
    --sintax_cutoff 0.5 \
    --tabbedout Path/to/output_file \
    --strand plus \
    --threads 1
```

### Assignation with usearch_global:

```bash
vsearch --usearch_global Path/to/fasta_file_after_step_2 \
    --db "$PATH_DIR_DB/$NAME_DB" \
    --id 0.7 \
    --iddef 2 \
    --strand plus \
    --threads 1 \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id+target \
    --maxaccepts 1 \
    --maxrejects 32 \
    --userout path/to/output_file
```




    
    