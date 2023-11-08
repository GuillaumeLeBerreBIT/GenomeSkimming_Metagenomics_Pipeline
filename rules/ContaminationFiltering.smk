#!/usr/bin/python3
############################# RULE #############################
# 
# All the reads gathered in one file can be filtered on contaminated reads. Contaminated reads filtered against the KRAKEN standard database.
# The KRAKEN standard database contains: Refeq archaea, bacteria, viral, plasmid, human1, UniVec_Core from latest update 6/5/2023.
# Can view the taxonomic classification reads in either a simplified report in *.txt or in the full *.kraken file which then each read is classified. 
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# Performing KRAKEN using Standard database to filter the reads on contamination
rule ContaminationFiltering:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_All_Reads_Concat.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )
    output:
        OutputReads = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/12_Contaminant_Kraken2/{SAMPLE}/ClassifiedContamination.kraken"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/kraken2.yaml"

    params:
        dbk2 = config["genome_skimming"]["KrakenCont"],
        smp_report = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/12_Contaminant_Kraken2/{SAMPLE}/ReportContamination.txt"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    shell:
        """
        kraken2 --use-names --threads 8 --db {params.dbk2} --report {params.smp_report} {input} > {output.OutputReads}
        """