############################# WORKFLOW #############################
# 
# The AllAssemblies_KrakenBuild_workflow will do the complete Genome Skimming pipeline. 
# 1) Assembling mitochondrial reference genomes, ribosomal repeats generating a Genbank file with the results summerized. 
# 2) Performing the Assembly free method to save all reads after contamination filtering of Species
# 3) Using the reads of each species saved in FASTA file >> To create a Custom database of Genome Skimming data with KRAKEN and BRACKEN Build    
#
###################################################################
#
############################# RULES - PROGRAMS #############################
include:
    '../rules/PreprocessingStrict.smk'
include:
    '../rules/MitoAssembly.smk'
include:
    '../rules/RibosomalAssembly.smk'
include:
    '../rules/PreAssemblies.smk'
include:
    '../rules/MitoAnnotation.smk'
include:
    '../rules/RibosomalAnnotation.smk'
include:
    '../rules/MergePE.smk'
include:
    '../rules/ContaminationFiltering.smk'
include:
    '../rules/ExtractKrakenReads.smk'
include:
    '../rules/BuildKraken2DB.smk'   
include: 
    '../rules/CreateGenbank.smk'    
include: 
    '../rules/CreateReferenceGenomeMapping.smk'

############################# RULE - RESULTING ANALYSIS #############################
# This functions as an "rule all" where the input is expected to be the finalized output. 
# By using touch can create an empty file to fake the output. 
rule Analysis:
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
                os.path.join(DATA_DIR_GS, "{PROJECT}/17_Reference_Mapping_MG/{TAXID}_{SPECIES}.fasta"),
            ],
            PROJECT = config["genome_skimming"]["project"],
            SAMPLE = config["genome_skimming"]["sample"],
            DB = config['genome_skimming']['KrakenDB'],
            TAXID = config['genome_skimming']['taxid'],
            SPECIES = config["genome_skimming"]["species"]
        )
    output:
        touch('AllAssembly_KrakenBuild_analysis.done')


