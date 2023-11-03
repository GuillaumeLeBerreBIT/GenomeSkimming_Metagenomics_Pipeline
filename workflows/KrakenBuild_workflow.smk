############################# WORKFLOW #############################
# 
# The KrakenBuild_workflow will with the KRAKEN Build create a custom database of all fasta files present after performing AssemblyFree method. 
# 1) Using the reads of each species saved in FASTA file >> To create a Custom database of Genome Skimming data with KRAKEN and BRACKEN Build 
#
###################################################################

############################# RULES - PROGRAMS #############################
include:
    '../rules/BuildKraken2DB.smk'

############################# RULE - RESULTING ANALYSIS #############################
# This functions as an "rule all" where the input is expected to be the finalized output. 
# By using touch can create an empty file to fake the output. 
rule Analysis:
    input:
        expand(
                [
                    
                os.path.join(DATA_DIR_GS, "{PROJECT}/15_Kraken_Databases/{DB}/database100mers.kraken"),
                    
                ],
                DATA_DIR = config["genome_skimming"]["datadir"],
                PROJECT = config["genome_skimming"]["project"],
                DB = config['genome_skimming']['KrakenDB']
            )
    output:
        touch('BuildKraken2DB_analysis.done')