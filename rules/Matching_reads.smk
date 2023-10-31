# #!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# 
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################
# 
rule ContaminationFiltering:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Unpaired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Unpaired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    output:
        kraken = expand(
            os.path.join(DATA_DIR_MG, "{project}/03_KrakenContaminant/{sample}_ClassifiedContamination.kraken"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
        ),
        unclassified = expand(
            os.path.join(DATA_DIR_MG, "{project}/03_KrakenContaminant/{sample}_Unclassified.fasta"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
        
    params:
        dbk2 = os.path.join(KRAKEN_STAND, "K2_Standard_DB"),
        smp_report = expand(
            os.path.join(DATA_DIR_MG, "{project}/03_KrakenContaminant/{sample}SimpleReport.txt"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )

    conda: 
        "../envs/kraken2.yaml"

    shell:
        """
        kraken2 --threads 16 --use-names --db {params.dbk2} --report {params.smp_report} --unclassified-out {output.unclassified} \
        {input.PAIRED_1} {input.UNPAIRED_1} {input.PAIRED_2} {input.UNPAIRED_2} > {output.kraken}
        """

rule ClassifyingReads:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/03_KrakenContaminant/{sample}_Unclassified.fasta"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/04_Classified_Kraken2/{sample}_Classified_Report.kraken"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    
    conda: 
        "../envs/kraken2.yaml"

    params:
        dbk2 = config["metagenomics"]["KrakenCustomDB"],
        kraken = expand(
            os.path.join(DATA_DIR_MG, "{project}/04_Classified_Kraken2/{sample}_Classified.kraken"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
            

    shell:
        """
        kraken2 --threads 16 --use-names --db {params.dbk2} --report {output} {input} > {params.kraken}
        """

