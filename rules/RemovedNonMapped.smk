#!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# 
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################

rule SamtoolsView:
    input:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/07_BWA_Mapped_Sequences/{sample}_Mapping.done"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/08_Sam_Filtered/{sample}_Samtools_view.done"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )

    conda:
        "../envs/samtools.yaml"

    params:
        InputMappings = expand(os.path.join(DATA_DIR_MG, "{project}/07_BWA_Mapped_Sequences/{sample}/"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        OutputFiltered = expand(os.path.join(DATA_DIR_MG, "{project}/08_Sam_Filtered/{sample}/"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        Sample = config["metagenomics"]["sample"]
    
    shell:
        """
        python3 scripts/RunSamtools.py -mf {params.InputMappings} -of {params.OutputFiltered}
        touch {output}
        """
        # Firstly ${sam##*/} will remove the file name from the path. 
        # {samfile%%.*} will split the name from the file extension
