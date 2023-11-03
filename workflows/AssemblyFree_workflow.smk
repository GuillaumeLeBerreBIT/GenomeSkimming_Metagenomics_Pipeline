############################# WORKFLOW #############################
# 
# The AssemblyFree_workflow will do AssemblyFree method of the Genome Skimming workflow. 
# 1) Performing the Assembly free method to save all reads after contamination filtering of Species
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
                os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/{TAXID}_{SAMPLE}.fasta")
            ],
            PROJECT = config["genome_skimming"]["project"],
            SAMPLE = config["genome_skimming"]["sample"],
            TAXID = config['genome_skimming']['taxid']
        )
    output:
        touch('AssemblyFree_analysis.done')

