#!/usr/bin/python3
############################# RULE #############################
# 
# This Snakemake rule is to perform the BRACKEN analysis from KRAKEN Classified reads. 
# Redistributing the reads for estimation of species abundances in metagenomic samples
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# 
rule Bracken:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/04_Classified_Kraken2/{SAMPLE}_Classified_Report.kraken"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )

    output:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/05_BRACKEN_Results/{SAMPLE}_Bracken_Classified.bracken"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    params:
        dbk2 = config["metagenomics"]["KrakenCustomDB"]

    conda: 
        "../envs/bracken.yaml"

    shell:
        """
        bracken -d {params.dbk2} -i {input} -o {output}
        """