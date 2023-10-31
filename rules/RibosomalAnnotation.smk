#!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# Annotation of the Ribosomak repeat regions. First ITSx is used to determine the 18S end, 5.8S start & end, 28S start region. 
# The annotated regions are togheter in a fasta file, this file is annotated again by Barrnap to determine the 
# start of 18S and end of 28S.   
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################
# 
rule ITSx:
    input:
        expand(os.path.join(DATA_DIR_GS, "{project}/05_GetOrganelle_Ribo_Results/{sample}/get_org.log.txt"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )
    output:
        expand(os.path.join(DATA_DIR_GS, "{project}/07_ITSx_Anno_Results/{sample}/ITSx_Anno.full.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )

    conda:
        "../envs/itsx.yaml"

    params:
        prefix = expand(os.path.join(DATA_DIR_GS, "{project}/07_ITSx_Anno_Results/{sample}/ITSx_Anno"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        ),
        # Use the * to take the file from GetOrganelle, which has a random output name --- > Only want the first file. 
        FTA = expand(os.path.join(DATA_DIR_GS, "{project}/05_GetOrganelle_Ribo_Results/{sample}/*1.1.*path_sequence.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )
        

    shell:
        """
        ITSx -t M -i {params.FTA} -o {params.prefix} --cpu 4 --save_regions all
        """

rule ConcatAllRegions:
    input:
        expand(os.path.join(DATA_DIR_GS, "{project}/07_ITSx_Anno_Results/{sample}/ITSx_Anno.full.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )
    output:
       expand(os.path.join(DATA_DIR_GS, "{project}/07_ITSx_Anno_Results/{sample}/{sample}_Concat_All_Regions_ITSx.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )

    conda:
        "../envs/itsx.yaml"

    params:
        prefix = expand(os.path.join(DATA_DIR_GS, "{project}/07_ITSx_Anno_Results/{sample}"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )
        
    shell:
        """
        cat {params.prefix}/*.SSU.fasta {params.prefix}/*.ITS1.fasta {params.prefix}/*.5_8S.fasta \
        {params.prefix}/*.ITS2.fasta {params.prefix}/*.LSU.fasta > {output}
        """

rule Barrnap:
    input:
        expand(os.path.join(DATA_DIR_GS, "{project}/07_ITSx_Anno_Results/{sample}/{sample}_Concat_All_Regions_ITSx.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )
    output:
        expand(os.path.join(DATA_DIR_GS, "{project}/08_Barrnap_Anno_Results/{sample}/{sample}_rDNA.gff"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )

    conda:
        "../envs/barrnap.yaml"

    params:
        FTA = expand(os.path.join(DATA_DIR_GS, "{project}/05_GetOrganelle_Ribo_Results/{sample}/*1.1.*path_sequence.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        ),
        RiboSeq = expand(os.path.join(DATA_DIR_GS, "{project}/08_Barrnap_Anno_Results/{sample}/{sample}_rDNA.fasta"),
        project = config['genome_skimming']['project'], sample = config['genome_skimming']['sample']
        )


    shell:
        """
        barrnap --thread 4 --outseq {params.RiboSeq} --kingdom euk {input} > {output}
        """