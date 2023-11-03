#!/usr/bin/python3
############################# RULE #############################
# 
# Using a python script to run BWA-MEM MEM & BWA-MEM INDEX, using One sample to iterate over all the
# files in '17_Reference_Mapping_MG'. 
# Commands run: 
# "bwa-mem2 index -p {file_index} {file_path}"
# "bwa-mem2 mem -t 8 {file_index} {args.PE_1} {args.PE_2} -o {args.MappingOutput}/{args.sample}_{file_splitted[0]}.sam"
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# Python script running BWA in the Shell
rule BWA_Index:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        #UNPAIRED_1 = expand(
        #    os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Unpaired.fastq"),
        #    PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
        #    ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        #UNPAIRED_2 = expand(
        #    os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Unpaired.fastq"),
        #    PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
        #    )        
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/07_BWA_Mapped_Sequences/{SAMPLE}_Mapping.done"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
        
    conda:
        "../envs/bwa-mem2.yaml"

    params:
        FolderReferenceFasta =  config["metagenomics"]["FolderReferenceFasta"],
        OutputFolderIndexes = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/06_Mapped_Indexes/"),
            PROJECT = config["metagenomics"]["project"]
            ),
        Sample = config["metagenomics"]["sample"],
        OutputMappings = expand(os.path.join(DATA_DIR_MG, "{PROJECT}/07_BWA_Mapped_Sequences/{SAMPLE}/"),
            PROJECT = config["metagenomics"]["project"], 
            SAMPLE = config["metagenomics"]["sample"]
            )
    
    threads: 8

    shell:
        """
        python3 scripts/RunBWA-MEM2.py -f {params.FolderReferenceFasta} -1 {input.PAIRED_1} -2 {input.PAIRED_2} -s {params.Sample} -io {params.OutputFolderIndexes} -mo {params.OutputMappings} -t {output}
        """
