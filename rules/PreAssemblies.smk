# #!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# With SPAdes will create all possible assemblies with only the paired files. 
# Filtering all preassmblies shorter then e.g. 200 bp.
# Needs to further optimized/worked on. 
#
###################################################################
#
############################# MODULES #############################
import os, glob

############################# RULES #############################
# The rule that invokes the SPAdes Assembly command to create all possible assemblies. 
rule SPAdesAssmebly:
    input:
        PAIRED_1 = os.path.join(DATA_DIR_GS, "{project}/02_Trimmomatic_Results/{sample}_for_paired.fastq"),
        UNPAIRED_1 = os.path.join(DATA_DIR_GS, "{project}/02_Trimmomatic_Results/{sample}_for_unpaired.fastq")
    output:
        GetOrg = os.path.join(DATA_DIR_GS, "{project}/06_SPAdes_PreA_Results/{sample}/scaffolds.fasta")

    conda:
        "../envs/spades.yaml"

    params:
        SPAdesFolder = os.path.join(DATA_DIR_GS, "{project}/06_SPAdes_PreA_Results/{sample}/")

    shell:
        """
        spades.py --isolate -1 {input.PAIRED_1} -2 {input.PAIRED_2} -o {params.SPAdesFolder}
        """


rule FindORFs:
    input:
    # Take the file with the resulting contigs.
        os.path.join(DATA_DIR_GS, "{project}/06_SPAdes_PreA_Results/{sample}/scaffolds.fasta")
    output:
        ""
    conda:
        ""
    shell:
        """
        
        """