#!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# This rule does the preprocessing of the FASTQ data. It does a Quality Control of the FASTQ data
# as well as Trimming of the reads removing Illumina adapters with Trimmomatic. 
# At last want to merge the PE-reads using NGmerge.
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# Trimmomatic removes the adapters from Illumina data. Reads matching in forward and reverse will be placed in a paired file. 
# unpaired reads will be placed in seperate files. 
rule Trimmomatic:
    input:
        FOR = expand(
            os.path.join(DATA_DIR_MG, "{SEQ_DIR}/{SAMPLE}{R1}"),
            SEQ_DIR = config["metagenomics"]["seqdata"], 
            SAMPLE = config["metagenomics"]["sample"],
            R1 = config["metagenomics"]["suffix_r1"]
            ),
        REV = expand(
            os.path.join(DATA_DIR_MG, "{SEQ_DIR}/{SAMPLE}{R2}"), 
            SEQ_DIR = config["metagenomics"]["seqdata"], 
            SAMPLE = config["metagenomics"]["sample"],
            R2 = config["metagenomics"]["suffix_r2"]
            )
    output:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Unpaired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Unpaired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )

    conda: 
        "../envs/trimmomatic.yaml"

    shell:
        """
        trimmomatic PE -threads 8 {input.FOR} {input.REV} {output.PAIRED_1} {output.UNPAIRED_1} {output.PAIRED_2} {output.UNPAIRED_2} \
        ILLUMINACLIP:adapters/NexteraPE-PE.fa:2:30:10 HEADCROP:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

# This is to check the quality of the reads using FastQC.
rule FastQ:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Unpaired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Unpaired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    
    output:
        PAIRED_1 = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Paired_fastqc.html"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_1 = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R1_Unpaired_fastqc.html"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Paired_fastqc.html"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_2 = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/{SAMPLE}_R2_Unpaired_fastqc.html"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )

    conda:
        "../envs/fastqc.yaml"

    params:
        FastqFol = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/02_FastQC_Results/"),
            PROJECT = config["metagenomics"]["project"]
            )

    shell:
        """
        fastqc --threads 4 {input.PAIRED_1} {input.UNPAIRED_1} {input.PAIRED_2} {input.UNPAIRED_2} -o {params.FastqFol}
        """
