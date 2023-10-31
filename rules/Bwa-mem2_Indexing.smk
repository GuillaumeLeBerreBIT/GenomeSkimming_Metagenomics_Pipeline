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


rule BWA_Index:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        #UNPAIRED_1 = expand(
        #    os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Unpaired.fastq"),
        #    project = config["metagenomics"]["project"], sample=config["metagenomics"]["sample"]
        #    ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        #UNPAIRED_2 = expand(
        #    os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Unpaired.fastq"),
        #    project = config["metagenomics"]["project"], sample=config["metagenomics"]["sample"]
        #    )        
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/07_BWA_Mapped_Sequences/{sample}_Mapping.done"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
        
    conda:
        "../envs/bwa-mem2.yaml"

    params:
        FolderReferenceFasta =  config["metagenomics"]["FolderReferenceFasta"],
        OutputFolderIndexes = expand(os.path.join(DATA_DIR_MG, "{project}/06_Mapped_Indexes/"),
            project = config["metagenomics"]["project"]
            ),
        Sample = config["metagenomics"]["sample"],
        OutputMappings = expand(os.path.join(DATA_DIR_MG, "{project}/07_BWA_Mapped_Sequences/{sample}/"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    
    threads: 8

    shell:
        """
        python3 scripts/RunBWA-MEM2.py -f {params.FolderReferenceFasta} -1 {input.PAIRED_1} -2 {input.PAIRED_2} -s {params.Sample} -io {params.OutputFolderIndexes} -mo {params.OutputMappings}
        touch {output}
        """
