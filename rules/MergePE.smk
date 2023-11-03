#!/usr/bin/python3
############################# RULE #############################
#  
# Want to merge the PE-reads using Fastq-join.
# Followed by concatenating all the reads Merged or not in one file to create one big DB.
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################

# Can merge the forward and reversed reads by using a tool called Fastq-join. The PE-reads are merged in a file.
# In the output have two more options to define the unpaired reads for R1 and R2. 
rule FastqJoin:
    input: 
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_for_paired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_back_paired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )
    
    output:
        MERGED = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_Merged.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        R1 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_R1.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        R2 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_R2.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/fastqjoin.yaml"

    shell:
        """
        fastq-join {input.PAIRED_1} {input.PAIRED_2} -o {output.R1} -o {output.R2} -o {output.MERGED}
        """

# Paste all the reads togheter in one file, merged and unmerged reads to have one FIle having all reads of the Species Genome Skim HT-Seq data. 
rule concatReads:
    input: 
        MERGED = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_Merged.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        R1 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_R1.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        R2 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_R2.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_for_unpaired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_back_unpaired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )
        
    
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_All_Reads_Concat.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    shell:
        """
        cat {input.MERGED} {input.R1} {input.R2} {input.UNPAIRED_1} {input.UNPAIRED_2} > {output}
        """