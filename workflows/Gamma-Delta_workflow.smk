############################# INTRODUCTION #############################
# 
# 
#
###################################################################

############################# RULES - PROGRAMS #############################
include:
    '../rules/PreprocessingStand.smk'
include:
    '../rules/Bwa-mem2_Indexing.smk'
include:
    '../rules/RemovedNonMapped.smk'
include:
    '../rules/Gamma-Delta_Algorithm.smk'
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
                os.path.join(DATA_DIR_MG, "{project}/08_Sam_Filtered/{sample}_Samtools_view.done")
                    
                ],
                project = config["metagenomics"]["project"],
                sample = config["metagenomics"]["sample"]                
            )
    output:
        touch('Gamma-Delta_Analysis.done')