# #!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# All the reads gathered in one file can be filtered on contaminated reads. Contaminated reads filtered against the KRAKEN standard database.
# The KRAKEN standard database contains: Refeq archaea, bacteria, viral, plasmid, human1, UniVec_Core from latest update 6/5/2023.
# Can view the taxonomic classification reads in either a simplified report in *.txt or in the full *.kraken file which then each read is classified. 
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################
# 
rule ContaminationFiltering:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_All_Reads_Concat.fastq"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{project}/12_Contaminant_Kraken2/{sample}/ClassifiedContamination.kraken"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/kraken2.yaml"

    params:
        dbk2 = "databases/K2_Standard_DB/",
        smp_report = expand(
            os.path.join(DATA_DIR_GS, "{project}/12_Contaminant_Kraken2/{sample}/ReportContamination.txt"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            )

    shell:
        """
        kraken2 --use-names --threads 8 --db {params.dbk2} --report {params.smp_report} {input} > {output}
        """