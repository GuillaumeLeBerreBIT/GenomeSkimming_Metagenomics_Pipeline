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
rule Bracken:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/04_Classified_Kraken2/{sample}_Classified_Report.kraken"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )

    output:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/05_BRACKEN_Results/{sample}_Bracken_Classified.bracken"), 
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    params:
        dbk2 = config["metagenomics"]["KrakenCustomDB"]

    conda: 
        "../envs/bracken.yaml"

    shell:
        """
        bracken -d {params.dbk2} -i {input} -o {output}
        """