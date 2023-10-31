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
rule CreateMappingFasta:
    input:
        MitosFile = expand(os.path.join(DATA_DIR_GS, "{project}/04_GetOrganelle_Mito_Results/{sample}/get_org.log.txt"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            ),
        RiboFile = expand(os.path.join(DATA_DIR_GS, "{project}/05_GetOrganelle_Ribo_Results/{sample}/get_org.log.txt"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            )


    output:
        expand(os.path.join(DATA_DIR_GS, "{project}/17_Reference_Mapping_MG/{taxid}_{species}.fasta"),
            project = config["genome_skimming"]["project"], 
            species = config["genome_skimming"]["species"],
            taxid = config["genome_skimming"]["taxid"]
            )

    conda:
        "../envs/biopython.yaml"

    params:
        MitosFolder = expand(os.path.join(DATA_DIR_GS, "{project}/04_GetOrganelle_Mito_Results/{sample}/"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            ),
        RiboFolder = expand(os.path.join(DATA_DIR_GS, "{project}/05_GetOrganelle_Ribo_Results/{sample}/"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            ),
        Species = config['genome_skimming']['species']

    shell:
        """
        python3 scripts/Combine_Mito_Ribo_Fasta.py -m {params.MitosFolder} -r {params.RiboFolder} -o {output} -s {params.Species}
        """