############################# WORKFLOW #############################
# 
# The AssemblyFree_KrakenBuild_workflow will do AssemblyFree method combined with the KRAKEN Build to create a 
# custom database of all fasta files present after performing AssemblyFree method. 
# 1) Performing the Assembly free method to save all reads after contamination filtering of Species
# 2) Using the reads of each species saved in FASTA file >> To create a Custom database of Genome Skimming data with KRAKEN and BRACKEN Build 
#
###################################################################
#
############################# RULES - PROGRAMS #############################
include:
    '../rules/PreprocessingStrict.smk'
include:
    '../rules/MergePE.smk'
include:
    '../rules/ContaminationFiltering.smk'
include:
    '../rules/ExtractKrakenReads.smk'
include:
    '../rules/BuildKraken2DB.smk'

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
                os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/database100mers.kraken"),
            ],
            PROJECT = config["genome_skimming"]["project"],
            SAMPLE = config["genome_skimming"]["sample"],
            DB = config['genome_skimming']['KrakenDB']
        )
    output:
        touch('AssemblyFree_KrakenBuild_analysis.done')