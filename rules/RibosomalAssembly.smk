# #!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# After trimming the reads and doing a FastQC can use the paired and unpaired reads for the ribosomal assembly.
# Using a customized seed and label database to assemble the ribosomal 18S-ITS1-5.8S-ITS2-28S region.    
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# 
rule GetOrganelleRibo:
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
        GetOrg = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/05_GetOrganelle_Ribo_Results/{SAMPLE}/get_org.log.txt"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/getorganelle.yaml"

    params:
        GetOrgFolder = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/05_GetOrganelle_Ribo_Results/{SAMPLE}/"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        Label_DB = "databases/GetOrg_Custom_DB/label_database.fasta",
        Seed_DB = "databases/GetOrg_Custom_DB/seed_database.fasta"

    threads: 8
    
    shell:
        """
        get_organelle_from_reads.py -1 {input.PAIRED_1} -2 {input.PAIRED_2} -u {input.UNPAIRED_1},{input.UNPAIRED_2} -t 8 \
        --reduce-reads-for-coverage inf --max-reads inf --genes {params.Label_DB} -s {params.Seed_DB} \
        -o {params.GetOrgFolder} -F anonym -R 10 --continue
        """