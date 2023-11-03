#!/usr/bin/python3
############################# RULE #############################
# 
# The trimmed reads can be used for the mitochondrial assembly performed by GetOrganelle. 
# 
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# This will perform the mitochondrial genome assembly
# The only file that has a consistent name is "get_org.log.txt" which will use as target output, since its the only file with a consistent name.
# Due to Snakemake creating directory added --continue otherwise it crashes here because it will create the folder itself.
# --overwrite completly overwrites the folder
rule GetOrganelleMito:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_for_paired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_for_unpaired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_back_paired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_back_unpaired.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )
    output:
        expand(os.path.join(DATA_DIR_GS, "{PROJECT}/04_GetOrganelle_Mito_Results/{SAMPLE}/get_org.log.txt"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/getorganelle.yaml"

    params:
        GetOrgFolder = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/04_GetOrganelle_Mito_Results/{SAMPLE}/"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    threads: 8
    
    shell:
        """
        get_organelle_config.py --add animal_mt
        get_organelle_from_reads.py -1 {input.PAIRED_1} -2 {input.PAIRED_2} -u {input.UNPAIRED_1},{input.UNPAIRED_2} -t 8 \
        --reduce-reads-for-coverage inf --max-reads inf -o {params.GetOrgFolder} -F animal_mt -R 10 --continue
        """
