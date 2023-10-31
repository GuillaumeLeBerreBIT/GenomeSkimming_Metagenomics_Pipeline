############################# INTRODUCTION #############################
# 
# 
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
                    
                os.path.join(DATA_DIR_GS, "{project}/15_Kraken_Databases/{DB}/database100mers.kraken"),
                    
                ],
                DATA_DIR = config["genome_skimming"]["datadir"],
                project = config["genome_skimming"]["project"],
                DB = config['genome_skimming']['KrakenDB']
            )
    output:
        touch('BuildKraken2DB_analysis.done')