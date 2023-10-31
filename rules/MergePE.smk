# #!/usr/bin/python3
############################# INTRODUCTION #############################
#  
# Want to merge the PE-reads using Fastq-join.
# Followed by concatenating all the reads Merged or not in one file to create one big DB.
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################

# Can merge the forward and reversed reads by using a tool called Fastq-join. The PE-reads are merged in a file.
# In the output have two more options to define the unpaired reads for R1 and R2. 
rule FastqJoin:
    input: 
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_GS, "{project}/01_Trimmomatic_Results/{sample}_for_paired.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_GS, "{project}/01_Trimmomatic_Results/{sample}_back_paired.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            )
    
    output:
        MERGED = expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_Merged.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        R1 = expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_R1.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        R2 = expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_R2.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
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
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_Merged.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        R1 = expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_R1.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        R2 = expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_R2.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_GS, "{project}/01_Trimmomatic_Results/{sample}_for_unpaired.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_GS, "{project}/01_Trimmomatic_Results/{sample}_back_unpaired.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            )
        
    
    output:
        os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_All_Reads_Concat.fastq")

    shell:
        """
        cat {input.MERGED} {input.R1} {input.R2} {input.UNPAIRED_1} {input.UNPAIRED_2} > {output}
        """