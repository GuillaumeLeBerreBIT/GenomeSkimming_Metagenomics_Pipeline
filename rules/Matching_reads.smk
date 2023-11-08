# #!/usr/bin/python3
############################# RULE #############################
# 
# Determining the taxonomic level for the reads from environmental samples. First will classify the reads on contamination database. 
# The unclassified reads will be used to match the reads against Kraken Database of Genome Skimming data
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# Do contamination filtering of the environmental samples. 
rule ContaminationFiltering:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R1_Unpaired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Paired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            ),
        UNPAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/01_Trimmomatic_Results/{SAMPLE}_R2_Unpaired.fastq"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    output:
        kraken = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/03_KrakenContaminant/{SAMPLE}_ClassifiedContamination.kraken"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
        ),
        unclassified = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/03_KrakenContaminant/{SAMPLE}_Unclassified.fasta"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
        
    params:
        dbk2 = config["metagenomics"]["KrakenCont"]
        smp_report = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/03_KrakenContaminant/{SAMPLE}SimpleReport.txt"),
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )

    conda: 
        "../envs/kraken2.yaml"

    shell:
        """
        kraken2 --threads 16 --use-names --db {params.dbk2} --report {params.smp_report} --unclassified-out {output.unclassified} \
        {input.PAIRED_1} {input.UNPAIRED_1} {input.PAIRED_2} {input.UNPAIRED_2} > {output.kraken}
        """

# Classify the reads against a custom created Genome Skim database. 
rule ClassifyingReads:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/03_KrakenContaminant/{SAMPLE}_Unclassified.fasta"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/04_Classified_Kraken2/{SAMPLE}_Classified_Report.kraken"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
    
    conda: 
        "../envs/kraken2.yaml"

    params:
        dbk2 = config["metagenomics"]["KrakenCustomDB"],
        kraken = expand(
            os.path.join(DATA_DIR_MG, "{PROJECT}/04_Classified_Kraken2/{SAMPLE}_Classified.kraken"), 
            PROJECT = config["metagenomics"]["project"], SAMPLE = config["metagenomics"]["sample"]
            )
            
    shell:
        """
        kraken2 --threads 16 --use-names --db {params.dbk2} --report {output} {input} > {params.kraken}
        """

