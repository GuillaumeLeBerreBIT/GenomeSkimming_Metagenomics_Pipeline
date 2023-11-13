#!/usr/bin/python3
############################# RULE #############################
# 
# Gamma-Delta Algorithm to determine the if a read maps on multiple species, see if the read can be
# classified to one species and if not not classify the read. 
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# Python script that runs the Gamma_Delta algorithm
rule GammaDelta:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/08_Sam_Filtered/{SAMPLE}_Samtools_view.done"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    output:
        summary = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/09_Gamma_Delta_Results/{SAMPLE}/{SAMPLE}_Gamma_Delta_Summary.csv"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )

    conda:
        "../envs/biopython.yaml"

    params:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        OutputFiltered = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/08_Sam_Filtered/{SAMPLE}/"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        assignment = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/09_Gamma_Delta_Results/{SAMPLE}/{SAMPLE}_Gamma_Delta_Reads_Assignment.csv"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        reads = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/09_Gamma_Delta_Results/{SAMPLE}/{SAMPLE}"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
        
    shell:
        """
        python3 scripts/gamma-delta.py -g 0.99 -d 0.98 -m {params.OutputFiltered} -r1 {params.PAIRED_1} -r2 {params.PAIRED_2} \
        -o {output.summary} -O {params.assignment} -F {params.reads}
        """