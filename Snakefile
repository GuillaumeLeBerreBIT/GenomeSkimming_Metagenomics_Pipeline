############################# INTRODUCTION #############################
# 
# The main Snakefile functions as possibility to choose the different Workflows for the analysis. 
# The Assembly method is to perfrom a miotchondrial, ribosomal and pre-assembly 
# There is also the possibility to perform a Assembly free method were the reads are matched
# based on k-mer index. 
#
###################################################################
#
############################# PARAMS #############################
configfile: "config.yaml"

DATA_DIR_GS = config["genome_skimming"]["datadir"]
DATA_DIR_MG = config["metagenomics"]["datadir"]

DB_DIR = config['mitos']['database']
DK_DIR = config["genome_skimming"]["dockerdir"]

# For the assembly free final name of output file -- > Multiple FASTA but with same TAXID ordered
TAXID = config['genome_skimming']['taxid']

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
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
                        os.path.join(DATA_DIR_GS, "{project}/08_Barrnap_Anno_Results/{sample}/{sample}_rDNA.gff"),
                        os.path.join(DATA_DIR_GS, "{project}/16_Genbank/{TAXID}_{sample}.done"),
                        os.path.join(DATA_DIR_GS, "{project}/17_Reference_Mapping_MG/{TAXID}_{species}.fasta"),
                    ],
                    project = config["genome_skimming"]["project"],
                    sample = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB'],
                    TAXID = config['genome_skimming']['taxid'],
                    species = config["genome_skimming"]["species"]
                )

    # Classify the reads using an Assembly free method which is based on k-mer matching of a certain length in a read.
    elif config["genome_skimming"]["workflow"] == "AssemblyFree":
        include:
            "workflows/AssemblyFree_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/14_Kraken2_Library_Reads/{TAXID}_{sample}.fasta"),
                        
                    ],
                    project = config["genome_skimming"]["project"],
                    sample = config["genome_skimming"]["sample"],
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
        
                        os.path.join(DATA_DIR_GS, "{project}/15_Kraken_Databases/{DB}/database100mers.kraken"),
                        
                    ],
                    project = config["genome_skimming"]["project"],
                    DB = config['genome_skimming']['KrakenDB']
                )

    # The full pipeline to start of FASTQ data from GenomeSkim HT-Seq data to use an AssemblFree method using the k-mer matching to create a GenomeSkim with all reads present in the GenomeSkim database (14_Kraken2_Library_Reads). 
    elif config["genome_skimming"]["workflow"] == "AssemblyFree_KrakenBuild":
        include:
            "workflows/AssemblyFree_KrakenBuild_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/15_Kraken_Databases/{DB}/database100mers.kraken")
                    ],
                    project = config["genome_skimming"]["project"],
                    sample = config["genome_skimming"]["sample"],
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
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
                        os.path.join(DATA_DIR_GS, "{project}/08_Barrnap_Anno_Results/{sample}/{sample}_rDNA.gff"),
                        os.path.join(DATA_DIR_GS, "{project}/14_Kraken2_Library_Reads/{TAXID}_{sample}.fasta"),
                        os.path.join(DATA_DIR_GS, "{project}/16_Genbank/{TAXID}_{sample}.done"),
                        os.path.join(DATA_DIR_GS, "{project}/17_Reference_Mapping_MG/{TAXID}_{species}.fasta"),
                    ],
                    project = config["genome_skimming"]["project"],
                    sample = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB'],
                    TAXID = config['genome_skimming']['taxid'],
                    species = config["genome_skimming"]["species"]
                )

    # This performs all the Assembly methods to either create a mtDNA, rDNA, Nuclar Assemblies as the AssemblyFree way of the k-mer matching of reads. Expanded to build the GenomeSkim database with it as well (Performing the complete GenomeSkim pipeline).
    elif config["genome_skimming"]["workflow"] == "AllAssemblies_KrakenBuild":
        include:
            "workflows/AllAssemblies_KrakenBuild_workflow.smk"
        rule all:
            input:
                expand(
                    [   
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
                        os.path.join(DATA_DIR_GS, "{project}/08_Barrnap_Anno_Results/{sample}/{sample}_rDNA.gff"),
                        os.path.join(DATA_DIR_GS, "{project}/16_Genbank/{TAXID}_{sample}.done"),
                        os.path.join(DATA_DIR_GS, "{project}/15_Kraken_Databases/{DB}/database100mers.kraken"),
                        os.path.join(DATA_DIR_GS, "{project}/17_Reference_Mapping_MG/{TAXID}_{species}.fasta")
                    ],
                    project = config["genome_skimming"]["project"],
                    sample = config["genome_skimming"]["sample"],
                    DB = config['genome_skimming']['KrakenDB'],
                    TAXID = config['genome_skimming']['taxid'],
                    species = config["genome_skimming"]["species"]
                )

elif config["pipeline"] == "Metagenomics":

    ## METAGENOMICS
    if config["metagenomics"]["workflow"] == "Kmer_Matching":
        include:
            "workflows/Kmer_matching_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R1_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R1_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R2_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R2_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/05_BRACKEN_Results/{sample}_Bracken_Classified.bracken")
                        
                    ],
                    project = config["metagenomics"]["project"],
                    sample = config["metagenomics"]["sample"]                
                )
    
    ## METAGENOMICS
    elif config["metagenomics"]["workflow"] == "Gamma_Delta":
        include:
            "workflows/Gamma-Delta_workflow.smk"
        rule all:
            input:
                expand(
                    [
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R1_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R1_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R2_Paired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/02_FastQC_Results/{sample}_R2_Unpaired_fastqc.html"),
                        os.path.join(DATA_DIR_MG, "{project}/08_Sam_Filtered/{sample}_Samtools_view.done")
                        
                    ],
                    project = config["metagenomics"]["project"],
                    sample = config["metagenomics"]["sample"]                
                )

else:
    raise Exception("Unknown workflow option: %s" % config["genome_skimming"]["workflow"])