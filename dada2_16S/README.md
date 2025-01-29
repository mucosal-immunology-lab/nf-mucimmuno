# DADA2 16S rRNA amplicon sequencing pre-processing

- [DADA2 16S rRNA amplicon sequencing pre-processing](#dada2-16s-rrna-amplicon-sequencing-pre-processing)
  - [Introduction](#introduction)
  - [Usage](#usage)
    - [Download the repository :open\_file\_folder:](#download-the-repository-open_file_folder)
    - [Create the conda environment :snake:](#create-the-conda-environment-snake)
    - [Folder structure](#folder-structure)
    - [Prepare your sample sheet :pencil:](#prepare-your-sample-sheet-pencil)
    - [Running the pipeline :running:](#running-the-pipeline-running)
      - [Customisation :gear:](#customisation-gear)
  - [Outputs](#outputs)

## Introduction

**nf-mucimmuno/dada2_16S** is a bioinformatics pipeline that can be used to run the popular DADA2 pre-processing pipeline for bacterial 16S rRNA amplicon sequencing data. It can handle multiple runs to generate a single unified output. It takes a samplesheet and either pre- or post-demultiplexed data (depending on what you have available), performs quality profiling, filtering and trimming, sequencing error modelling, sample inference, and merging of paired ends. From there, it combines all runs together, removes chimeras, assigns SILVA taxonomy, and generate a *de novo* phylogenetic tree using RAxML.

<div align="center">
<img src="../assets/images/nf-mucimmuno_dada2_16S.png" width=80%>
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
SUBFOLDER="dada2_16S"

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

Then run the script to retrieve the repository into a new folder called `dada2_16S`, which will house your workflow files and results.

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
mamba activate nextflow-dada2
```

### Folder structure

Because you specify the full directory path for your raw input data, you can technically house them however and wherever you like. However, below is an example of how to store runs whether you have demultiplexed files or the raw sequencing files.

- `run01` requires demultiplexing, and as such has a single file for each of `R1.fastq.gz`, `R2.fastq.gz`, and `Index.fastq.gz`. It also importantly contains a the barcode mapping file which **must** include `barcode_to_sample` in its name.
- `run02` however is already demultiplexed, and therefore only requires the forward (`R1`) and reverse (`R2`) reads.
  - Critically, the pipeline looks for an existing `demultiplexed` folder to identify whether it can skip demultiplexing.

```bash
dada2_16S/
    └── data/
        ├── run01/
        │   ├── barcode_to_sample_run01.txt
        │   ├── Index.fastq.gz
        │   ├── R1.fastq.gz
        │   └── R2.fastq.gz
        └── run02/
            └── demultiplexed/
                ├── sample1-R1.fastq.gz
                ├── sample1-R2.fastq.gz
                ├── sample2-R1.fastq.gz
                ├── sample2-R2.fastq.gz
                └── ...
```

### Prepare your sample sheet :pencil:

This pipeline requires a sample sheet to identify where your sequencing data is located. You can also change the name of your run to something some specific if you desire from the original folder name.

Your sample sheet should look as follows, **ensuring you use the exact column names as below**.

```bash
run_name,folder_path
run_01,/path_to_pipeline/dada2_16S/data/run01
run_02,/path_to_pipeline/dada2_16S/data/run02
```

An example is provided [here](./sample_sheet.csv).

### Running the pipeline :running:

Now you can run the pipeline. You will need to set up a parent job to run each of the individual jobs &ndash; this can be either an interactive session, or an sbatch job. For example:

```bash
# Start an interactive session with minimal resources
smux n --time=3-00:00:00 --mem=16GB --ntasks=1 --cpuspertask=2 -J nf-STARsolo
```

Make sure you alter the `nextflow.config` file to provide the path to your sample sheet, unless it is `./samplesheet.csv` which is the default for the cluster profile. Stay within the top `cluster` profile section to alter settings for Slurm-submitted jobs.

Inside your interactive session, be sure to activate your `nextflow-scrnaseq` environment from above. Then, **inside the scRNAseq folder**, begin the pipeline using the following command (ensuring you use the `cluster` profile to make use of the Slurm workflow manager).

```bash
# Activate conda environment
mamba activate nextflow-dada2

# Begin running the pipeline
nextflow run dada2_pipeline.nf -resume -profile cluster
```

#### Customisation :gear:

There are some customisation options that are available within the `nextflow.config` file. While the defaults should be suitable for those with access to the M3 MASSIVE cluster genomics partition, for those without access, of for those who require different amounts of resources, there are ways to change these.

To adjust the `cluster` profile settings, stay within the appropriate section at the top of the file.

***Parameters***

| Option | Description |
| --- | --- |
| sample_sheet | The file path to your sample sheet |
| outdir | A new folder name to be created for your results |
| *filter_and_trim*.truncLen |  |
| *filter_and_trim*.maxEE |  |
| *filter_and_trim*.trimLeft |  |
| *filter_and_trim*.truncQ |  |
| *filter_and_trim*.maxN |  |
| *assign_taxonomy*.trainSet_link |  |
| *assign_taxonomy*.trainSet_file |  |
| *assign_taxonomy*.assignSpecies_link |  |
| *assign_taxonomy*.assignSpecies_file |  |

The DADA2-ready SILVA taxonomic data is sourced from [here](https://zenodo.org/records/14169026). If new updates become available, feel free to update these values (or let us know so we can update).

***Process***

| Option | Description |
| --- | --- |
| executor | The workload manager (default: `'slurm'`) |
| conda | The conda environment to use (default: `'./environment.yaml'`) |
| queueSize | The maximum number of jobs to be submitted at any time (default: `12`) |
| submitRateLimit | The rate allowed for job submission &ndash; either a number of jobs per second (e.g. 20sec) or a number of jobs per time period (e.g. 20/5min) (default: `'1/2sec'`) |
| memory | The maximum global memory allowed for Nextflow to use (default: `'320 GB'`) |
| *DEMULTIPLEX*.memory | Memory for DEMULTIPLEX step to use (default: `'52 GB'`) |
| *DEMULTIPLEX*.cpus | Number of CPUs for DEMULTIPLEX step to use (default: `'8'`) |
| *DEMULTIPLEX*.clusterOptions | Specific cluster options for DEMULTIPLEX step to use (default: `'--time=4:00:00 --partition=genomics --partition=genomics'`) |
| *ERROR_MODEL*.memory | Memory for ERROR_MODEL step to use (default: `'52 GB'`) |
| *ERROR_MODEL*.cpus | Number of CPUs for ERROR_MODEL step to use (default: `'8'`) |
| *ERROR_MODEL*.clusterOptions | Specific cluster options for ERROR_MODEL step to use (default: `'--time=4:00:00 --partition=genomics --partition=genomics'`) |
| *MERGE_PAIRED_ENDS*.memory | Memory for MERGE_PAIRED_ENDS step to use (default: `'52 GB'`) |
| *MERGE_PAIRED_ENDS*.cpus | Number of CPUs for MERGE_PAIRED_ENDS step to use (default: `'8'`) |
| *MERGE_PAIRED_ENDS*.clusterOptions | Specific cluster options for MERGE_PAIRED_ENDS step to use (default: `'--time=4:00:00 --partition=genomics --partition=genomics'`) |
| *REMOVE_CHIMERAS*.memory | Memory for REMOVE_CHIMERAS step to use (default: `'52 GB'`) |
| *REMOVE_CHIMERAS*.cpus | Number of CPUs for REMOVE_CHIMERAS step to use (default: `'8'`) |
| *REMOVE_CHIMERAS*.clusterOptions | Specific cluster options for REMOVE_CHIMERAS step to use (default: `'--time=4:00:00 --partition=genomics --partition=genomics'`) |
| *ASSIGN_TAXONOMY*.memory | Memory for ASSIGN_TAXONOMY step to use (default: `'52 GB'`) |
| *ASSIGN_TAXONOMY*.cpus | Number of CPUs for ASSIGN_TAXONOMY step to use (default: `'8'`) |
| *ASSIGN_TAXONOMY*.clusterOptions | Specific cluster options for ASSIGN_TAXONOMY step to use (default: `'--time=4:00:00 --partition=genomics --partition=genomics'`) |
| *DE_NOVO_PHYLO_TREE*.memory | Memory for v step to use (default: `'160 GB'`) |
| *DE_NOVO_PHYLO_TREE*.cpus | Number of CPUs for DE_NOVO_PHYLO_TREE step to use (default: `'24'`) |
| *DE_NOVO_PHYLO_TREE*.clusterOptions | Specific cluster options for DE_NOVO_PHYLO_TREE step to use (default: `'--time=24:00:00'`) |

The *de novo* phylogenetic tree step has been given extra resources and time by default, as this can take a very long time to run. You can always experiment with reducing these resources if you want; you can just rerun this step afterwards with more resources if it doesn't work.

## Outputs

Several outputs will be copied from their respective Nextflow `work` directories to the output folder of your choice (default: `results`).

The main outputs of interest for downstream processing are the inputs to create your `phyloseq` object.

| Output | Description |
| --- | --- |
| `seqtab_nochim.rds` | The clean counts table matrix |
| `taxonomy_species.rds` | The SILVA-assigned taxonomy table |
| `tree.rds` | The *de novo* phylogenetic tree |

There are also a collection of quality control and filtering plots available in this directory too. Below is an example of the output structure for running a single run.

```bash
dada2_16S/
    ├── run01/
    │   └── filtering_report/
    │       └── run01_filtered/
    │           ├── error_model.pdf
    │           ├── error_model.pdf.rds
    │           ├── error_model.rds
    │           ├── filtering_report.pdf
    │           ├── filtering_report.pdf.rds
    │           ├── sample1-R1.fastq
    │           ├── sample1-R2.fastq
    │           ├── sample2-R1.fastq
    │           ├── sample2-R2.fastq
    │           └── ...
    ├── aggregated_error_plots.pdf
    ├── aggregated_filtering_plots.pdf
    ├── aggregated_quality_profiles.pdf
    ├── length_distribution.pdf
    ├── length_distribution.rds
    ├── read_tracking_report.pdf
    ├── seqtab_nochim.rds
    ├── seqtab_nochim.txt
    ├── seqtab.rds
    ├── taxonomy_species.rds
    ├── taxonomy_species.txt
    ├── taxonomy.rds
    ├── taxonomy.txt
    └── tree.rds
```