# GenomeSkimming_Metagenomics_Pipeline 

This Snakemake workflow offers a comprehensive solution for the analysis of high-throughput sequencing data obtained through genome skimming. It enables the application of various assembly methods and metagenomic analyses, providing flexibility and robustness to researchers in the field.

Assembly Methods
Our workflow provides two distinct assembly methods:

### Assembly Method:

Utilizing GetOrganelle, this method constructs a reference database comprising complete mitochondrial genomes and ribosomal repeats. The assembly process ensures the accurate representation of genetic material for downstream analysis.

### Assembly-Free Method:

In this approach, the workflow employs KRAKEN to generate a k-mer database from all sequencing reads, post-contamination filtering. This KRAKEN database serves as a fundamental resource for subsequent metagenomic analysis.

Metagenomic Analysis
The workflow supports two distinct modes of metagenomic analysis:

### Reference-Based Analysis:

Metagenomic sequencing reads are aligned to the self-created reference mitochondrial genomes and ribosomal repeats. The alignment results are subsequently employed in a taxonomic classification process using the Gamma-Delta algorithm. This method allows for a detailed exploration of the metagenomic composition.

### K-mer-Based Analysis:

Alternatively, researchers have the option to perform k-mer matching of the reads using KRAKEN, followed by BRACKEN. This approach enables the efficient estimation of taxonomic abundance within the sample, especially suited for diverse and complex metagenomic datasets.

This README provides a detailed guide to the setup, configuration, and execution of the workflow, allowing users to harness the full potential of their genome skimming data for assembly and metagenomic insights.

## Installation

Before you can start using the workflow need to have Conda, Singularity and Snakemake installed to perform the analysis. 

### Conda

Snakemake requires to have a conda installation on your system. The preferred conda distribution is mambaforge since it has the required python commands & mamba which is a very fast installation. 
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
When the download of miniconda is complete can start the local installation. The installation can be completly up to you where to install it on your system. 
```
bash Miniconda3-latest-Linux-x86_64.sh
```
Snakemake will make use of the conda environment to install envs provided by yaml files. The configuration of the .condarc file is very important for Snakemake!
```
conda config --add channels conda-forge
conda config --add channels bioconda
```
The resulting .condarc file looks like this. 
```
auto_activate_base: false
channels:
  - bioconda
  - conda-forge
  - defaults
```
WATCH OUT! When having problems with installing conda environments through snakemake, check if the channel_priority is set to strict. Can delete it from the .condarc file and try again. 

### Snakemake

Snakemake can be installed using a conda. Mamba is possible to use as well if installed.
```
conda create -c conda-forge -c bioconda -n snakemake snakemake=7.25.0
```
Can check if snakemake is installed. The version should be 7.25 if installed using previous command. 
```
conda activate snakemake
snakemake --help
snakemake --version
```

### Singularity

When cloning the whole directory, can find in resources the `singularity-ce.tar.gz` file. Can extract the directory in the main folder. 
```
tar xvzf resources/singularity-ce.tar.gz -C .
```
After extracted need to set the path to the singularity-ce folder so snakemake can use the installation to perform singularity installations. 
```
nano ~/.bashrc
export PATH=$PATH:/path/to/singularity-ce/bin

source ~/.bashrc
```

## Usage

The configuration settings govern the execution of workflows, dictating which processes are initiated and where data is stored or utilized. Users can specify folder names for data and results, databases, workflows, sample names, and more.

Within the Snakemake workflow, the primary distinction lies between two pipelines: `Genome_Skimming` for creating reference genomes, and `Metagenomics` for taxonomic classification of reads.
```
pipeline: "Genome_Skimming" 
```

The `Genome_Skimming` pipeline offers two options for creating a reference database: assembling reference genomes or employing a k-mer-based approach. The sample name for Genome Skimming is derived from the filename, split into two parts - the common name between 'R1' and 'R2' read names (`sample`), and the suffixes of the forward and reverse reads (`suffix_r1` and `suffix_r2`).

The `datadir` configuration allows users to define the path to the folder containing all sequencing data and metadata.

For workflows involving the `Assembly` analysis and the use of singularity, it is mandatory to mount an empty directory to the project folder. Ensure that the specified directory for mounting exists in the folder where the Snakefile is executed. Command line arguments can be used to specify the mounted path to the `project` directory, saving processed metadata within it.

The `KrakenCont` configuration allows users to set the path to the Kraken database used for contaminant filtering during the AssemblyFree analyses. Users can either utilize a custom Kraken database or download one from [Pre Build Databases](https://benlangmead.github.io/aws-indexes/k2). Once downloaded can decompress the database and define the path in the configuration file. 

The `workflow` configuration determines the type of analysis under the Genome Skimming workflow. Options include

- Assembly: To perform the Assemblies of mtDNA, rDNA and Nuclear DNA. 
- AssemblyFree: Perfrom the k-mer matching of the Genome Skim reads using KRAKEN2
- KrakenBuildDB: Build a custom KRAKEN2 database from all gathered reads after Aseembly Free
- AssemblyFree_KrakenBuild: The complete pipeline to perfrom k-mer matching of the reads, followed by the creation of a custom database. 
- AllAssemblies: Perfrom the Assembly methods for mtDNA, rDNA and Nuclear DNA togheter with the k-mer matching og the GenomeSkim reads.
- AllAssemblies_KrakenBuild: Perfrom the Assembly methods for mtDNA, rDNA and Nuclear DNA togheter with the k-mer matching og the GenomeSkim reads. Finalized with the creation of the GenomeSkim database with KRAKEN2. 

The `species` configuration sets the name of the animal sample, while `taxid` configures the Taxonomy Identification linked to the species name using NCBI Taxonomy identifiers.

During the `KrakenBuildDB` workflow, a custom Kraken database is built, and its name can be defined under `KrakenDB`.
```
genome_skimming:
  sample: 'MB_GS_O_ophiura_S36'

  suffix_r1: '_R1_001.fastq.gz'
  suffix_r2: '_R2_001.fastq.gz'
  
  datadir: "/home/genomics/gleberre/01_Research_BAR_ZAND/02_ZAND/" 

  seqdata: "00_ZAND_RAW_SEQ/01_SEQ_GS"

  project: "01_ZAND_GS"

  dockerdir: "/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/mnt"

  KrakenCont: "/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/databases/K2_Standard_DB/"

  workflow: "AssemblyFree" 
  
  species: "Ophuria_ophuria"

  taxid: "72673" 

  KrakenDB: "MacrobenthosDB"
```

Under the `mitos` configuration, users can specify the path to the MITOS reference database (database) and the genetic codon table (`gencode`). WATCH OUT! The used database, being the `refseq89m` is not added to the path. Mitos requires 2 arguments one being the path to the used reference database and one for the name of the used reference database. Other database are possible from either [MITOS](https://zenodo.org/records/2683856) or [MITOS2](https://zenodo.org/records/4284483). 

Under `gencode` can define the genetic codon table want to use. 

- 5 == Invertebrate
- 2 == Vertebrate
- 4 == Mold,
- 14 == Alternative flatworm
- 9 == Echinoderm
- 13 == Ascidia

```
mitos:
  database: "/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/databases/MITOS_DB/"
  gencode: 5
```

The `genbank` configurations only need to be set/changed when running Assembly, AllAssemblies, AllAssemblies_KrakenBuild analysis since it will be used to create a Genbank file. The LOCUS name can be defined for the mitochondrial genome `locus_mt` or ribosomal DNA repeats `locus_rdna`. Can set the bioproject name for mitochondrial genbank `bioproject_mt` and ribosomal genbank `bioproject_rdna`. There are other optionts that can be set such as `sample_type`, `taxonomy` where the name has underscore for spaces, `organism_header` can have spaces in the name. Can set the `mol_type` which by default is DNA. The species name defined in the SOURCE of the genbank is configured at `organisms`. 
```
genbank:

  locus_mt: "Notomastus_latericeus_mtDNA"
  locus_rdna: "Notomastus_latericeus_rDNA"

  bioproject_mt: "Mitochondrial_Assembly"
  bioproject_rdna: "Ribosomal_Assembly"

  sample_type: "ILVO - RMG"
  taxonomy: "Notomastus_latericeus"
  organism_header: "Notomastus latericeus"
  mol_type: "DNA"

  organism: "Notomastus_latericeus"
``` 

Under `metagenomics`, users can define configurations for the metagenomics analysis, similar to those for Genome Skimming.

The `workflow` is what will define what analysis to you will run under the Metagenomics workflow. 
  
- Kmer_Matching: To perform k-mer matching of the reads against a "Custom" created database
- Gamma_Delta: Perform the mapping of reads against reference mitochondrial organelle and ribosomal repeats
- Gamma_Delta_Kmer: Will do both the mapping as well as the matching of the reads

Under `KrakenCustomDB` can define the path to the (Custom) Kraken database want to use to classify the metagenomic sequencing reads. 

The `FolderReferenceFasta` is the path to the folder `17_Reference_Mapping_MG` containing all the Reference (multi)FASTA files which are used for mapping reads using BWA-MEM under the `Gamma-Delta` analysis.  

```
metagenomics:
  
  sample: 'BAR_MG_Field_3_S28'

  suffix_r1: '_R1_001.fastq.gz'
  suffix_r2: '_R2_001.fastq.gz'
  
  workflow: "Gamma_Delta"
  
  seqdata: "00_BAR_RAW_SEQ/02_SEQ_MG"

  project: "02_BAR_MG"
  
  datadir: "/home/genomics/gleberre/01_Research_BAR_ZAND/01_BAR/" 

  KrakenStand: "/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/databases/K2_Standard_DB/"

  KrakenCustomDB: "/home/genomics/gleberre/01_Research_BAR_ZAND/01_BAR/01_BAR_GS/15_Kraken_Databases/GenomeSkimDB"

  FolderReferenceFasta: "/home/genomics/gleberre/01_Research_BAR_ZAND/01_BAR/01_BAR_GS/17_Reference_Mapping_MG"
```

### Commands 

To run the workflow can use the following command which works when using any workflow defined under Genome Skimming and metagenomics. Running the argument `--use-conda` makes sure the conda environments are used. Because MITOS is a docker image need `--use-singularity` to run A paired with `--singularity-args` followed by the command argument want to use when executing singularity. The command option used is `-B` which is to specify the `/path/to/directory/to/mnt:/path/to/container/`. At last the number of cores are added which if using `KrakenBuild` workflows require at least 16 cores. 
```
snakemake --use-conda --use-singularity --singularity-args "-B /home/genomics/gleberre/01_Research_BAR_ZAND/01_BAR/01_BAR_GS:/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/mnt" --cores 24
```

For workflows that do not involve singularity/docker images, a shorter command using only conda environments suffices.
```
snakemake --use-conda --cores 16 
```