#!/usr/bin/python3
############################# RULE #############################
# 
# Running a python script that iterates over the mapped results of a sample. 
# Filtering for Unmapped, Not primary alignment and Supplementary alignment
# Command:
# "samtools view -@ 8 -S -F 2308 {sam_file_path} -o {args.OutputFilteredSam}/{basename[0]}_Filtered.sam"
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
#A python script that uses Samtools in the Shell
rule SamtoolsView:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/07_BWA_Mapped_Sequences/{SAMPLE}_Mapping.done"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/08_Sam_Filtered/{SAMPLE}_Samtools_view.done"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )

    conda:
        "../envs/samtools.yaml"

    params:
        InputMappings = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/07_BWA_Mapped_Sequences/{SAMPLE}/"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        OutputFiltered = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/08_Sam_Filtered/{SAMPLE}/"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        Sample = config["metagenomics"]["sample"]
    
    shell:
        """
        python3 scripts/RunSamtools.py -mf {params.InputMappings} -of {params.OutputFiltered}
        touch {output}
        """
        # Firstly ${sam##*/} will remove the file name from the path. 
        # {samfile%%.*} will split the name from the file extension
