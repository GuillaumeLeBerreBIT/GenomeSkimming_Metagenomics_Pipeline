############################# INTRODUCTION #############################
# 
# The main Snakefile functions as possibility to choose the different Workflows for the analysis high throughput sequencing reads. 
# The Assembly method is to assemble a mitochondrial reference genome and ribosomal repeats of high-throughput sequencing data.  
# There are also workflows that perform an Assembly-Free method, reads are filtered against a Database of Bacteria, Refseq Archea, Viral, ...
# the resulting reads that were unclassified are saved in a Fasta file 'TAXID_SPECIES.fasta'. Those files can be used to craete a custom KRAKEN database
# 
# There are also workflows to do metagenomics analysis.
# Can match reads against the Kraken custom Genome Skimming database to classifiy reads from environmental samples, combinad with BRACKEN to 
# re-distributing reads in the taxonomic tree. 
# Reads can also be classified using created reference genomes through assembly methods, where the read will be mapped against.
# Based on how well the read maps on different reference species it will classify the read to one specific species. 
#
###################################################################
#
############################# PARAMS #############################
configfile: "config.yaml"
# Paths to the directory containing Raw sequencing data, Project folder (All generated metadata)
DATA_DIR_GS = config["genome_skimming"]["datadir"]
DATA_DIR_MG = config["metagenomics"]["datadir"]
# Database used for MITOS Annotation
DB_DIR_MITOS = config['mitos']['database']
# Path to the docker directory used for the MITOS Annotation
# Used to mount a directory to a folder in folder where the main Snakemake file is present
DK_DIR = config["genome_skimming"]["dockerdir"]

# The Taxonomy ID is added to the output files such as Genbank, Reference file for mapping and 
# after contamination filtering the TAXID is added to the Kraken reads that were unclassified in 
# the contamination filtering as well the output file names of Kraken  classification
TAXID = config['genome_skimming']['taxid']
# Path to the KRAKEN2 standard database used for contamination filtering
KRAKEN_STAND = config["metagenomics"]["KrakenStand"]

############################# WORKFLOWS #############################
if config["pipeline"] == "Genome_Skimming":
    # This workflow will perfrom an Assembly of the Mitochondrial genome, Ribosomal and Nuclear free assembly. 
    if config["genome_skimming"]["workflow"] == "Assembly":
        include:
            "workflows/Assembly_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/06_MITOS_Results/{SAMPLE}/AnnotationMitosDone.txt"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/08_Barrnap_Anno_Results/{SAMPLE}/{SAMPLE}_rDNA.gff"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/16_Genbank/{TAXID}_{SAMPLE}.done"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/17_Reference_Mapping_MG/{TAXID}_{SPECIES}.fasta"),
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB'],
                    TAXID = config['genome_skimming']['taxid'],
                    SPECIES = config["genome_skimming"]["species"]
                )

    # Classify the reads using an Assembly free method which is based on k-mer matching of a certain length in a read.
    elif config["genome_skimming"]["workflow"] == "AssemblyFree":
        include:
            "workflows/AssemblyFree_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/{TAXID}_{SAMPLE}.fasta"),
                        
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"],
                    TAXID = config['genome_skimming']['taxid']
                )

    # Using the reads after the contamination from the Kraken2 standard database, to create the GenomeSkim database of unclassified reads. 
    elif config["genome_skimming"]["workflow"] == "KrakenBuildDB":
        include:
            "workflows/KrakenBuild_workflow.smk"
        rule all:
            input:
                expand(
                    [
        
                        os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/database100mers.kraken"),
                        
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    DB = config['genome_skimming']['KrakenDB']
                )

    # The full pipeline to start of FASTQ data from GenomeSkim HT-Seq data to use an AssemblFree method using the k-mer matching 
    # to create a Custom Genome Skim database (14_Kraken2_Library_Reads). 
    elif config["genome_skimming"]["workflow"] == "AssemblyFree_KrakenBuild":
        include:
            "workflows/AssemblyFree_KrakenBuild_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/database100mers.kraken")
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB']
                )

    # This performs all the Assembly methods to either create a mtDNA, rDNA, Nuclar Assemblies as the AssemblyFree way of the k-mer matching of reads. 
    elif config["genome_skimming"]["workflow"] == "AllAssemblies":
        include:
            "workflows/AllAssemblies_workflow.smk"
        rule all:
            input:
                expand(
                    [   
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/06_MITOS_Results/{SAMPLE}/AnnotationMitosDone.txt"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/08_Barrnap_Anno_Results/{SAMPLE}/{SAMPLE}_rDNA.gff"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/{TAXID}_{SAMPLE}.fasta"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/16_Genbank/{TAXID}_{SAMPLE}.done"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/17_Reference_Mapping_MG/{TAXID}_{SPECIES}.fasta"),
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB'],
                    TAXID = config['genome_skimming']['taxid'],
                    SPECIES = config["genome_skimming"]["species"]
                )

    # This performs all the Assembly methods to either create a mtDNA, rDNA, (Nuclar Assemblies) as the AssemblyFree way of the k-mer matching of reads. 
    # Expanded to build the GenomeSkim database with KRAKEN as well (Performing the complete GenomeSkim pipeline).
    elif config["genome_skimming"]["workflow"] == "AllAssemblies_KrakenBuild":
        include:
            "workflows/AllAssemblies_KrakenBuild_workflow.smk"
        rule all:
            input:
                expand(
                    [   
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/06_MITOS_Results/{SAMPLE}/AnnotationMitosDone.txt"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/08_Barrnap_Anno_Results/{SAMPLE}/{SAMPLE}_rDNA.gff"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/16_Genbank/{TAXID}_{SAMPLE}.done"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/database100mers.kraken"),
                        os.path.join(DATA_DIR_GS, "{PROJECT}/17_Reference_Mapping_MG/{TAXID}_{SPECIES}.fasta")
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB'],
                    TAXID = config['genome_skimming']['taxid'],
                    SPECIES = config["genome_skimming"]["species"]
                )
                
## METAGENOMICS
elif config["pipeline"] == "Metagenomics":

    # This performs the k-mer matching of reads against a custom created KRAKEN Genome Skimming database. 
    if config["metagenomics"]["workflow"] == "Kmer_Matching":
        include:
            "workflows/Kmer_matching_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/05_BRACKEN_Results/{SAMPLE}_Bracken_Classified.bracken")
                        
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"]                
                )
    
    # This performs the mapping of reads against a self-created custom reference mitochondrial genome and ribosomal repeat, 
    # determining the species identification based on a Gamma-Delta algorithm to assign reads to species level.
    elif config["metagenomics"]["workflow"] == "Gamma_Delta":
        include:
            "workflows/Gamma-Delta_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/09_Gamma_Delta_Results/{SAMPLE}_Gamma_Delta.csv")
                        
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"]                
                )
    # This performs the k-mer matching of reads against a custom created KRAKEN Genome Skimming database. Combined with 
    # the mapping of reads against a self-created custom reference mitochondrial genome and ribosomal repeat
    elif config["metagenomics"]["workflow"] == "Gamma_Delta_Kmer":
        include:
            "workflows/Gamma-Delta_Kmer_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/09_Gamma_Delta_Results/{SAMPLE}_Gamma_Delta.csv"),
                        os.path.join(DATA_DIR_MG, "{PROJECT}/05_BRACKEN_Results/{SAMPLE}_Bracken_Classified.bracken")
                        
                    ],
                    PROJECT = config["genome_skimming"]["project"],
                    SAMPLE = config["genome_skimming"]["sample"]                
                )

else:
    # If non of the defining pipelines is chosen. 
    raise Exception("Unknown workflow option: %s" % config["genome_skimming"]["workflow"])