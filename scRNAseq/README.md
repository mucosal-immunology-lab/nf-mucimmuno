# Pre-processing of scRNAseq FASTQ files


- [Pre-processing of scRNAseq FASTQ files](#pre-processing-of-scrnaseq-fastq-files)
  - [Introduction](#introduction)
  - [Usage](#usage)
    - [Download the repository :open\_file\_folder:](#download-the-repository-open_file_folder)
    - [Create the conda environment :snake:](#create-the-conda-environment-snake)
    - [Prepare the genome :dna:](#prepare-the-genome-dna)
      - [Human genome files :man::woman:](#human-genome-files-manwoman)
      - [Mouse genome files :mouse:](#mouse-genome-files-mouse)
    - [Prepare your sample sheet :pencil:](#prepare-your-sample-sheet-pencil)
    - [Running the pipeline :running:](#running-the-pipeline-running)
      - [Customisation :gear:](#customisation-gear)
  - [Outputs](#outputs)
    - [Collected export files :package:](#collected-export-files-package)
    - [Reports :page\_facing\_up:](#reports-page_facing_up)
    - [STARsolo :star:](#starsolo-star)


## Introduction

**nf-mucimmuno/scRNAseq** is a bioinformatics pipeline that can be used to run quality control steps and alignment to a host genome using STARsolo. It takes a samplesheet and FASTQ files as input, performs FastQC, trimming and alignment, and produces an output `.tar.gz` archive containing the collected outputs from STARsolo, ready for further processing downstream in R. MultiQC is run on the FastQC outputs both before and after TrimGalore! for visual inspection of sample quality &ndash; output `.html` files are collected in the results.

<div align="center">
<img src="../assets/images/nf-mucimmuno_scRNAseq-01.png" width=80%>
</div>

## Usage

### Download the repository :open_file_folder:

This repository contains the relevant Nextflow workflow components, including a conda environment and submodules, to run the pipeline. To retrieve this repository alone, run the `retrieve_me.sh` script above.

:warning: Git `sparse-checkout` is required to retrieve just the **nf-mucimmuno/scRNAseq** pipeline. It was only introduced to Git in version 2.27.0, so ensure that the loaded version is high enough (or that there is a version loaded on the cluster at all). As of July 2024, the M3 MASSIVE cluster has version 2.38.1 available. :warning:

```bash
# Check git version
git --version

# Load git module if not loaded or insufficient version
module load git/2.38.1
```

First, create a new bash script file.

```bash
# Create and edit a new file with nano
nano retrieve_me.sh
```

Add the contents to the file, save, and close.

```bash
#!/bin/bash

# Define variables
REPO_URL="https://github.com/mucosal-immunology-lab/nf-mucimmuno"
REPO_DIR="nf-mucimmuno"
SUBFOLDER="scRNAseq"

# Clone the repository with sparse checkout
git clone --no-checkout $REPO_URL
cd $REPO_DIR

# Initialize sparse-checkout and set the desired subfolder
git sparse-checkout init --cone
git sparse-checkout set $SUBFOLDER

# Checkout the files in the subfolder
git checkout main

# Move the folder into the main folder and delete the parent
mv $SUBFOLDER ../
cd ..
rm -rf $REPO_DIR

echo "Subfolder '$SUBFOLDER' has been downloaded successfully."
```

Then run the script to retrieve the repository into a new folder called `scRNAseq`, which will house your workflow files and results.

```bash
# Run the script
bash retrieve_me.sh
```

### Create the conda environment :snake:

To create the conda environment, use the provided environment `.yaml` file. Then activate it to access required functions.

```bash
# Create the environment
mamba env create -f environment.yaml

# Activate the environment
mamba activate nextflow-scrnaseq
```

### Prepare the genome :dna:

Create a new folder somewhere to store your genome files. Enter the new folder, and run the relevant code depending on your host organism. Run these steps in an interactive session with ~48GB RAM and 16 cores, or submit them as an sbatch job.

:warning: Please check if these are already available somewhere before regenerating them yourself! :warning:

STAR should be loaded already via the conda environment for the genome indexing step. We will set `--sjdbOverhang` to 78 to be suitable for use with the longer `R2` FASTQ data resulting from BD Rhapsody single cell sequencing. This may require alteration for other platforms.

#### Human genome files :man::woman:

```bash
#!/bin/bash
VERSION=111
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz
gunzip *
```

Then use STAR to prepare the genome index.

```bash
#!/bin/bash
VERSION=111
STAR \
    --runThreadN 16 \
    --genomeDir "STARgenomeIndex78/" \
    --runMode genomeGenerate \
    --genomeFastaFiles "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" \
    --sjdbGTFfile "Homo_sapiens.GRCh38.$VERSION.gtf" \
    --sjdbOverhang 78
```

#### Mouse genome files :mouse:

```bash
#!/bin/bash
VERSION=111
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/mus_musculus/Mus_musculus.GRCm39.$VERSION.gtf.gz
gunzip *
```

Then use STAR to prepare the genome index.

```bash
#!/bin/bash
VERSION=111
STAR \
    --runThreadN 16 \
    --genomeDir "STARgenomeIndex78/" \
    --runMode genomeGenerate \
    --genomeFastaFiles "Mus_musculus.GRCm39.dna_sm.primary_assembly.fa" \
    --sjdbGTFfile "Mus_musculus.GRCm39.$VERSION.gtf" \
    --sjdbOverhang 78
```

### Prepare your sample sheet :pencil:

This pipeline requires a sample sheet to identify where your FASTQ files are located, and which cell label sequences (CLS) are being utilised.

More information about the CLS tags used with BD Rhapsody single-cell RNAseq library preparation can be found here:

* [BD Rhapsody Sequence Analysis Pipeline &ndash; User's Guide](https://www.bdbiosciences.com/content/dam/bdb/marketing-documents/products-pdf-folder/software-informatics/rhapsody-sequence-analysis-pipeline/Rhapsody-Sequence-Analysis-Pipeline-UG.pdf)
* [BD Rhapsody Cell Label Structure &ndash; Python Script](https://bd-rhapsody-public.s3.amazonaws.com/CellLabel/rhapsody_cell_label.py.txt)

The benefit of providing the name of the CLS bead versions in the sample sheet is that you can combine runs that utilise different beads together in the same workflow. Keep in mind that if you do this though, there may be some bead-related batch effects to address and correct downstream &ndash; it is always important to check for these effects when combining sequencing runs in any case. The current options are:

| CLS option | Description |
|---|---|
| BD_Original | The original BD rhapsody beads and linker sequences |
| BD_Enhanced_V1 | First version of enhanced beads with polyT and 5prime capture oligo types, shorter linker sequences, longer polyT, and 0-3 diversity insert bases at the beginning of the sequence |
| BD_Enhanced_V2 | Same structure as the enhanced (V1) beads, but with increased CLS diversity (384 vs. 96) |

Your sample sheet should look as follows, **ensuring you use the exact column names as below**. Remember that on the M3 MASSIVE cluster, you need to use the **full file path** &ndash; relative file paths don't usually work.

```bash
sample,fastq_1,fastq_2,CLS
CONTROL_S1,CONTROL_S1_R1.fastq.gz,CONTROL_S1_R2.fastq.gz,BD_Enhanced_V2
CONTROL_S2,CONTROL_S2_R1.fastq.gz,CONTROL_S1_R2.fastq.gz,BD_Enhanced_V2
TREATMENT_S1,TREATMENT_S1_R1.fastq.gz,TREATMENT_S1_R2.fastq.gz,BD_Enhanced_V2
```

An example is provided in [`data/samplesheet_test`](./data/samplesheet_test.csv).

### Running the pipeline :running:

Now you can run the pipeline. You will need to set up a parent job to run each of the individual jobs &ndash; this can be either an interactive session, or an sbatch job. For example:

```bash
# Start an interactive session with minimal resources
smux n --time=3-00:00:00 --mem=16GB --ntasks=1 --cpuspertask=2 -J nf-STARsolo
```

Make sure you alter the `nextflow.config` file to provide the paths to both your sample sheet and the prepared genome index folder. Stay within the top `cluster` profile section to alter setting for Slurm-submitted jobs.

Inside your interactive session, be sure to activate your `nextflow-scrnaseq` environment from above. Then, **inside the scRNAseq folder**, begin the pipeline using the following command (ensuring you use the `cluster` profile to make use of the Slurm workflow manager).

```bash
# Begin running the pipeline
nextflow run process_raw_reads.nf -resume -profile cluster
```

#### Customisation :gear:

There are several customisation options that are available within the `nextflow.config` file. While the defaults should be suitable for those with access to the M3 MASSIVE cluster genomics partition, for those without access, of for those who require different amounts of resources, there are ways to change these.

To adjust the `cluster` profile settings, stay within the appropriate section at the top of the file.

***Parameters***

Visit [STAR documentation](https://github.com/alexdobin/STAR/) for explanations of all available options for STARsolo.

| Option | Description |
|---|---|
| samples_csv | The file path to your sample sheet |
| outdir | A new folder name to be created for your results |
| *trimgalore*.quality | The minimum quality before a sequence is truncated (default: `20`) |
| *trimgalore*.length | The minimum length of a sequence in order to be retained (default: `43`) |
| *trimgalore*.adapter | A custom adapter sequence for the R1 sequences (default: `'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'`) |
| *trimgalore*.adapter2 | A custom adapter sequence for the R2 sequences (default: `'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'`) |
| *starsolo*.index | The file path to the prepared genome index folder |
| *starsolo*.soloType | The type of single-cell RNAseq (default: `'CB_UMI_Complex'`) |
| *starsolo*.soloCBmatchWLtype | The method for matching cell barcodes to the whitelist files (default: `'EditDist_2'`) |
| *starsolo*.soloUMIdedup | The type of UMI deduplication (default: `'1MM_CR'`) |
| *starsolo*.soloUMIfiltering | The type of UMI filtering for reads uniquely mapping to genes (default: `'MultiGeneUMI_CR'`) |
| *starsolo*.soloCellFilter | The method type and parameters for cell filtering (default: `'EmptyDrops_CR'`) |
| *starsolo*.soloMultiMappers | The counting method for reads mapping for multiple genes (default: `'EM'`) |

***Process***

These settings relate to resource allocation and cluster settings. FASTQC and TRIMGALORE steps can take longer than 4 hours for typical single-cell RNAseq file, and therefore the default option is to run these steps on the `comp` partition.

| Option | Description |
| --- | --- |
| executor | The workload manager (default: `'slurm'`) |
| conda | The conda environment to use (default: `'./environment.yaml'`) |
| queueSize | The maximum number of jobs to be submitted at any time (default: `12`) |
| submitRateLimit | The rate allowed for job submission &ndash; either a number of jobs per second (e.g. 20sec) or a number of jobs per time period (e.g. 20/5min) (default: `'1/2sec'`) |
| memory | The maximum global memory allowed for Nextflow to use (default: `'320 GB'`) |
| *FASTQC*.memory | Memory for FASTQC step to use (default: `'80 GB'`) |
| *FASTQC*.cpus | Number of CPUs for FASTQC step to use (default: `8`) |
| *FASTQC*.clusterOptions | Specific cluster options for FASTQC step (default: `'--time=8:00:00'`) |
| *TRIMGALORE*.memory | Memory for TRIMGALORE step to use (default: `'80 GB'`) |
| *TRIMGALORE*.cpus | Number of CPUs for TRIMGALORE step to use (default: `8`) |
| *TRIMGALORE*.clusterOptions | Specific cluster options for TRIMGALORE step (default : `'--time=8:00:00'`) |
| *STARSOLO*.memory | Memory for STARSOLO step to use (default: `'80 GB'`) |
| *STARSOLO*.cpus | Number of CPUs for STARSOLO step to use (default: `12`) |
| *STARSOLO*.clusterOptions | Specific cluster options for STARSOLO step (default : `'--time=4:00:00 --partition=genomics --qos=genomics'`) |
| *COLLECT_EXPORT_FILES*.memory | Memory for COLLECT_EXPORT_FILES step to use (default: `'32 GB'`) |
| *COLLECT_EXPORT_FILES*.cpus | Number of CPUs for COLLECT_EXPORT_FILES step to use (default: `8`) |
| *COLLECT_EXPORT_FILES*.clusterOptions | Specific cluster options for COLLECT_EXPORT_FILES step (default : `'--time=4:00:00 --partition=genomics --qos=genomics'`) |

## Outputs

Several outputs will be copied from their respective Nextflow `work` directories to the output folder of your choice (default: `results`).

### Collected export files :package:

The main output will be a single archive file called `export_files.tar.gz` that you will take for further downstream pre-processing. It contains STARsolo outputs for each sample, with the respective subfolders described below.

### Reports :page_facing_up:

Within the `reports` folder, you will find the MultiQC outputs from pre- and post-trimming.

### STARsolo :star:

Contains the outputs for each sample from STARsolo, including various log files and package version information.

The main output of interest here is a folder called `{sample}.Solo.out`, which houses subfolders called `Gene`, `GeneFull_Ex50pAS`, and `Velocyto`. It is this main folder for each sample that is added to `export_files.tar.gz`.
* As you will use the gene count data from `GeneFull_Ex50pAS` downstream, it is a good idea to check the `Summary.csv` within this folder for each sample to ensure mapping was successful.
  * One of the key values to inspect is `Reads With Valid Barcodes`, which should be >0.8 (indicating at least 80% of reads had valid barcodes).
  * If you note that this value is closer to 0.02 (i.e. ~2% had valid barcodes), you should double-check to make sure you specified the correct BD Rhapsody beads version. For instance, if you specified `BD_Enhanced_V1` but actually required `BD_Enhanced_V2`, the majority of your reads will not match the whitelist, and therefore the reads will be considered invalid.

**Folder structure**

Below is an example of the output structure for running one sample. The STARsolo folder would contain additional samples as required.

```bash
scRNAseq
└── results/
    ├── export_files.tar.gz
    ├── reports/
    │   ├── pretrim_multiqc_report.html
    │   └── posttrim_multiqc_report.html
    └── STARsolo/
        └── sample1/
            ├── sample1.Solo.out/
            │   ├── Gene/
            │   │   ├── filtered/
            │   │   │   ├── barcodes.tsv.gz
            │   │   │   ├── features.tsv.gz
            │   │   │   └── matrix.mtx.gz
            │   │   ├── raw/
            │   │   │   ├── barcodes.tsv.gz
            │   │   │   ├── features.tsv.gz
            │   │   │   ├── matrix.mtx.gz
            │   │   │   └── UniqueAndMult-EM.mtx.gz
            │   │   ├── Features.stats
            │   │   ├── Summary.csv
            │   │   └── UMIperCellSorted.txt
            │   ├── GeneFull_Ex50pAS/
            │   │   ├── filtered/
            │   │   │   ├── barcodes.tsv.gz
            │   │   │   ├── features.tsv.gz
            │   │   │   └── matrix.mtx.gz
            │   │   ├── raw/
            │   │   │   ├── barcodes.tsv.gz
            │   │   │   ├── features.tsv.gz
            │   │   │   ├── matrix.mtx.gz
            │   │   │   └── UniqueAndMult-EM.mtx.gz
            │   │   ├── Features.stats
            │   │   ├── Summary.csv
            │   │   └── UMIperCellSorted.txt
            │   ├── Velocyto/
            │   │   ├── filtered/
            │   │   │   ├── ambiguous.mtx.gz
            │   │   │   ├── barcodes.tsv.gz
            │   │   │   ├── features.tsv.gz
            │   │   │   ├── spliced.mtx.gz
            │   │   │   └── unspliced.mtx.gz
            │   │   ├── raw/
            │   │   │   ├── ambiguous.mtx.gz
            │   │   │   ├── barcodes.tsv.gz
            │   │   │   ├── features.tsv.gz
            │   │   │   ├── spliced.mtx.gz
            │   │   │   └── unspliced.mtx.gz
            │   │   ├── Features.stats
            │   │   └── Summary.csv
            │   └── Barcodes.stats
            ├── sample1.Log.final.out
            ├── sample1.Log.out
            ├── sample1.Log.progress.out
            └── versions.yml
```