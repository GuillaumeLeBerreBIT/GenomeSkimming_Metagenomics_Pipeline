############################# INTRODUCTION #############################
# 
# This file will contain all the rules necessary to perform the analysis of the Assembly Free method. 
# The rules are split up in different files for efficient combining in different pipelines. 
#
###################################################################

############################# RULES - PROGRAMS #############################
include:
    '../rules/PreprocessingStrict.smk'
include:
    '../rules/MergePE.smk'
include:
    '../rules/ContaminationFiltering.smk'
include:
    '../rules/ExtractKrakenReads.smk'

############################# RULE - RESULTING ANALYSIS #############################
# This functions as an "rule all" where the input is expected to be the finalized output. 
# By using touch can create an empty file to fake the output. 
rule Analysis:
    input:
        expand(
            [
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/14_Kraken2_Library_Reads/{TAXID}_{sample}.fasta")
            ],
            project = config["genome_skimming"]["project"],
            sample = config["genome_skimming"]["sample"],
            TAXID = config['genome_skimming']['taxid']
        )
    output:
        touch('AssemblyFree_analysis.done')

