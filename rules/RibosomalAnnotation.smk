#!/usr/bin/python3
############################# RULE #############################
# 
# Annotation of the Ribosomal repeat regions. First ITSx is used to determine the 18S end, 5.8S start & end, 28S start region. 
# The annotated regions are concatted togheter in a fasta file, this file is annotated again by Barrnap to determine the 
# start of 18S and end of 28S.   
#
###################################################################
# 
############################# MODULES #############################
import os

############################# RULES #############################
# 
rule ITSx:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/05_GetOrganelle_Ribo_Results/{SAMPLE}/get_org.log.txt"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )
            
    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/07_ITSx_Anno_Results/{SAMPLE}/ITSx_Anno.full.fasta"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )

    conda:
        "../envs/itsx.yaml"

    params:
        prefix = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/07_ITSx_Anno_Results/{SAMPLE}/ITSx_Anno"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            ),
        # Use the * to take the file from GetOrganelle, which has a random output name --- > Only want the first file. 
        FTA = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/05_GetOrganelle_Ribo_Results/{SAMPLE}/*1.1.*path_sequence.fasta"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )
        

    shell:
        """
        ITSx -t M -i {params.FTA} -o {params.prefix} --cpu 4 --save_regions all
        """

rule ConcatAllRegions:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/07_ITSx_Anno_Results/{SAMPLE}/ITSx_Anno.full.fasta"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )

    output:
       expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/07_ITSx_Anno_Results/{SAMPLE}/{SAMPLE}_Concat_All_Regions_ITSx.fasta"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )

    conda:
        "../envs/itsx.yaml"

    params:
        prefix = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/07_ITSx_Anno_Results/{SAMPLE}"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )
        
    shell:
        """
        cat {params.prefix}/*.SSU.fasta {params.prefix}/*.ITS1.fasta {params.prefix}/*.5_8S.fasta \
        {params.prefix}/*.ITS2.fasta {params.prefix}/*.LSU.fasta > {output}
        """

rule Barrnap:
    input:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/07_ITSx_Anno_Results/{SAMPLE}/{SAMPLE}_Concat_All_Regions_ITSx.fasta"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )

    output:
        expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/08_Barrnap_Anno_Results/{SAMPLE}/{SAMPLE}_rDNA.gff"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )

    conda:
        "../envs/barrnap.yaml"

    params:
        RiboSeq = expand(
            os.path.join(DATA_DIR_GS, "{PROJECT}/08_Barrnap_Anno_Results/{SAMPLE}/{SAMPLE}_rDNA.fasta"),
            PROJECT = config['genome_skimming']['project'], SAMPLE = config['genome_skimming']['sample']
            )


    shell:
        """
        barrnap --thread 4 --outseq {params.RiboSeq} --kingdom euk {input} > {output}
        """