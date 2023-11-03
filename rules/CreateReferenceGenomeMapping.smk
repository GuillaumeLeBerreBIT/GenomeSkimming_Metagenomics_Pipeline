# #!/usr/bin/python3
############################# RULE #############################
# 
# Will combine the reference sequences from mitochondrial genome as well as the ribosomal repeat in one file
# The file can then be used for the mapping (index) the environmental samples. 
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################
# Combine FASTA sequences from the Assembly
rule CreateMappingFasta:
    input:
        MitosFile = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/04_GetOrganelle_Mito_Results/{SAMPLE}/get_org.log.txt"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        RiboFile = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/05_GetOrganelle_Ribo_Results/{SAMPLE}/get_org.log.txt"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            )

    output:
        expand(os.path.join(DATA_DIR_GS, "{PROJECT}/17_Reference_Mapping_MG/{TAXID}_{SPECIES}.fasta"),
            PROJECT = config["genome_skimming"]["project"], 
            SPECIES = config["genome_skimming"]["species"],
            TAXID = config["genome_skimming"]["taxid"]
            )

    conda:
        "../envs/biopython.yaml"

    params:
        MitosFolder = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/04_GetOrganelle_Mito_Results/{SAMPLE}/"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        RiboFolder = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/05_GetOrganelle_Ribo_Results/{SAMPLE}/"),
            PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
            ),
        Species = config['genome_skimming']['species']

    shell:
        """
        python3 scripts/Combine_Mito_Ribo_Fasta.py -m {params.MitosFolder} -r {params.RiboFolder} -o {output} -s {params.Species}
        """