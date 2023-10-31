############################# INTRODUCTION #############################
# 
# 
#
###################################################################

############################# RULES - PROGRAMS #############################
include:
    '../rules/PreprocessingStand.smk'
include:
    '../rules/Matching_reads.smk'
include: 
    '../rules/Bracken.smk'
############################# RULE - RESULTING ANALYSIS #############################
# 
rule Analysis:
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
    output:
        touch('Kmer_matching_Analysis.done')
