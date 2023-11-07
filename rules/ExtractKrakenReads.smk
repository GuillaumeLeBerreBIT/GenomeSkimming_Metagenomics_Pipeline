#!/usr/bin/python3
############################# RULE #############################
# 
# Can extract the reads that were unclassified == non contaminated reads. 
# Using a script provided by one of the contributors of KRAKEN to extract the reads on certain taxonomic and more. 
# At last for creating a custom DB need to add the taxonomic ID to the headers.
# It needs to match the |kraken:taxid|XXXX
# 
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################

# Add the Taxonomic Identifier to the headers of the reads from the GenomeSkim
rule Extract_Kraken2_Reads:
    input:
        OutputReads = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/12_Contaminant_Kraken2/{SAMPLE}/ClassifiedContamination.kraken"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        AllReads = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/03_Fastq_join_Results/{SAMPLE}/{SAMPLE}_All_Reads_Concat.fastq"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/13_Extract_Kraken2_Reads/{SAMPLE}/{SAMPLE}_Unclassified.fasta"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/biopython.yaml"

    shell:
        """
        python3 scripts/extract_kraken_reads.py -k {input.OutputReads} -s {input.AllReads} -t 0 -o {output}
        """

# Add the Taxonomic Identifier to the headers of the reads from the GenomeSkim
rule AddTaxonomicID:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/13_Extract_Kraken2_Reads/{SAMPLE}/{SAMPLE}_Unclassified.fasta"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/{TAXID}_{SAMPLE}.fasta"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"], TAXID = config['genome_skimming']['taxid']
            )

    conda:
        "../envs/biopython.yaml"

    params:
        taxid = config['genome_skimming']['taxid']

    shell:
        """
        python3 scripts/Add_Taxid.py -t {params.taxid} {input} {output}
        """