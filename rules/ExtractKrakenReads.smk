# #!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# Can extract the reads that were unclassified == non contaminated reads. 
# Using a script provided by one of the contributors of KRAKEN to extract the reads on certain taxonomic and more. 
# At last for creating a custom DB need to add the taxonomic ID to the headers.
# It needs to match the |kraken:taxid|XXXX
# 
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################
# Take out the unclassified reads that will be used to create the Genome Skim database. 
rule ExtractKrakenReads:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{project}/12_Contaminant_Kraken2/{sample}/ClassifiedContamination.kraken"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{project}/13_Extract_Kraken2_Reads/{sample}/Extracted_{sample}.fasta"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            )

    conda:
        "../envs/biopython.yaml"

    params:
        MERGED = expand(
            os.path.join(DATA_DIR_GS, "{project}/03_Fastq_join_Results/{sample}/{sample}_All_Reads_Concat.fastq"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            )

    shell:
        """
        python3 scripts/extract_kraken_reads.py -k {input} -s1 {params.MERGED} -t 0 -o {output}
        """

# Add the Taxonomic Identifier to the headers of the reads from the GenomeSkim
rule AddTaxonomicID:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{project}/13_Extract_Kraken2_Reads/{sample}/Extracted_{sample}.fasta"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{project}/14_Kraken2_Library_Reads/{TAXID}_{sample}.fasta"),
            project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"], TAXID = config['genome_skimming']['taxid']
            )

    conda:
        "../envs/biopython.yaml"

    params:
        taxid = config['genome_skimming']['taxid']

    shell:
        """
        python3 scripts/Add_Taxid.py -t {params.taxid}  {input} {output}
        """