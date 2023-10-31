#!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# 
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################

rule GammaDelta:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/08_Sam_Filtered/{sample}_Samtools_view.done"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/09_Gamma_Delta_Results/{sample}_Gamma_Delta.csv"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )

    params:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        OutputFiltered = expand(os.path.join(DATA_DIR_MG, "{project}/08_Sam_Filtered/{sample}/"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
        
    shell:
        """
        python3 scripts/gamma-delta.py -g 0.99 -d 0.98 -m {params.OutputFiltered} -r1 {params.PAIRED_1} -r2 {params.PAIRED_2} -o {output}
        """