#!/usr/bin/python3
############################# RULE #############################
# 
# At last can create a own database using the unclassified reads after the filtering of the contamination. 
# Firstly need to add all the Fasta files created for the GenomeSkim DB to the library which can be done with the kraken2-build --add-library.
# Even tho multiple fasta of the same species are present, it is dependant on the taxonomic classifier added in the ExtractKrakenReads.smk -- > AddTaxonomicID.
# After adding all the FASTA files can download the taxonomic classification with kraken2-build --download-taxonomy 
# When previous steps are completed can start building the custom KRAKEN DB with kraken2-build --build
# At last can clean up the database to remove unnecessary files.
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# Adding all the FASTA files to the library for the KRAKEN DB creation. 
rule AddLibrary:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/{TAXID}_{SAMPLE}.fasta"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"], TAXID = config['genome_skimming']['taxid']
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/BuildLibrary_{DB}.done"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )

    conda:
        "../envs/kraken2.yaml"

    params:
        KrakenFolder = expand(
                            os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/"),
                            PROJECT = config["genome_skimming"]["project"]
                            ),
        KrakenPathDB = expand(
                            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/"),
                            PROJECT = config["genome_skimming"]["project"]
                            ),
        DBname = config['genome_skimming']['KrakenDB']

    threads: 16

    shell:
        """
        for file in {params.KrakenFolder}/*.fasta
        do
            kraken2-build --threads 16 --add-to-library $file --db {params.KrakenPathDB}/{params.DBname}
        done
        touch {output}
        """

# Download the taxonomic data for the reads present in the library of the KRAKEN DB. 
rule AddTaxid:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/BuildLibrary_{DB}.done"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )

    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/AddTaxonomy_{DB}.done"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )

    conda:
        "../envs/kraken2.yaml"

    params:
        KrakenPathDB = expand(
                        os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/"),
                        PROJECT = config["genome_skimming"]["project"]
                        ),
        DBname = config['genome_skimming']['KrakenDB']
    
    threads: 16

    shell:
        """
        kraken2-build --threads 16 --use-ftp --download-taxonomy --db {params.KrakenPathDB}/{params.DBname}
        touch {output}
        """

# Create the KRAKEN database by building using all the taxonomic and library data. 
# Remove metadata that is not needed anymore. 
rule BuildKrakenDB:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/AddTaxonomy_{DB}.done"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/taxo.k2d"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )
        
    conda:
        "../envs/kraken2.yaml"

    params:
        KrakenPathDB = expand(
                        os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/"),
                        PROJECT = config["genome_skimming"]["project"]
                        ),
        DBname = config['genome_skimming']['KrakenDB']

    threads: 16

    # Will skip the clean step since need the library files for BRACKEN
    # kraken2-build --threads 16 --clean --db {params.KrakenPathDB}/{params.DBname}
    shell:
        """
        kraken2-build --threads 16 --build --db {params.KrakenPathDB}/{params.DBname}
        """

# After building the KRAKEN database can directly build the BRACKEN database since the 'database100mers.kmer_distrib' is needed to do BRACKEN
rule BrackenBuild:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/taxo.k2d"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/database100mers.kraken"),
            DB = config['genome_skimming']['KrakenDB'], PROJECT = config["genome_skimming"]["project"]
            )
    conda:
        "../envs/bracken.yaml"

    params:
        KrakenPathDB = expand(
                        os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/"),
                        PROJECT = config["genome_skimming"]["project"]
                        ),
        DBname = config['genome_skimming']['KrakenDB']
    
    threads: 16

    shell:
        """
        bracken-build -d {params.KrakenPathDB}/{params.DBname} -t 16
        """
    

