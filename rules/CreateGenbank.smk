#!/usr/bin/python3
############################# RULE #############################
# 
# When the Assembly for the mitochondrial genome and the annotation is generated. Can use the resulting FASTA and GFF file
# for creating the Genbank file. All kind of options are defined in the config which are parsed on the command line to create the wanted Genbank file.
# 
###################################################################
#
############################# MODULES #############################
import os

############################# RULES #############################
# Using the GFF and FASTA file to convert -- > Genbank file
rule CreateGenbank:
    input:
        AnnoMitos = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/06_MITOS_Results/{SAMPLE}/AnnotationMitosDone.txt"),
        PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
        ),
        AssemblyGetOrg = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/04_GetOrganelle_Mito_Results/{SAMPLE}/get_org.log.txt"),
        PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
        )

    output:
        expand(os.path.join(DATA_DIR_GS, "{PROJECT}/16_Genbank/{TAXID}_{SAMPLE}.done"),
        PROJECT = config["genome_skimming"]["project"], 
        SAMPLE = config["genome_skimming"]["sample"],
        TAXID = config['genome_skimming']['taxid']
        )

    conda:
        "../envs/biopython.yaml"

    params:
        FastaFile = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/04_GetOrganelle_Mito_Results/{SAMPLE}/*1.1.path_sequence.fasta"),
        PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
        ),
        GffFile = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/06_MITOS_Results/{SAMPLE}/"),
        PROJECT = config["genome_skimming"]["project"], SAMPLE = config["genome_skimming"]["sample"]
        ),
        GenbankFile = expand(os.path.join(DATA_DIR_GS, "{PROJECT}/16_Genbank/{TAXID}_{SAMPLE}.gb"),
        PROJECT = config["genome_skimming"]["project"], 
        SAMPLE = config["genome_skimming"]["sample"],
        TAXID = config['genome_skimming']['taxid']
        ),
        
        ## HEADER
        locus = config['genbank']['locus'],
        bioproject = config['genbank']['bioproject'],
        sample_type = config['genbank']['sample_type'],
        taxonomy = config['genbank']['taxonomy'],
        organism_header = config['genbank']['organism_header'],
        mol_type = config['genbank']['mol_type'],
        topology = config['genbank']['topology'],
        
        ## FEATURE TABLE
        # SOURCE
        organism = config['genbank']['organism'],
        organelle = config['genbank']['organelle']
    
    # /usr/bin/bash: -c: line 3: syntax error near unexpected token `(' >> Need to provided "", especially when spaces involved passed through the config file. 
    shell:
        """
        first_matching_file=$(ls -1 {params.FastaFile} | head -n 1)
        python3 scripts/Fasta_GFF_To_GBK.py -f $first_matching_file -g "{params.GffFile}" -gbk "{params.GenbankFile}" \
        -oi "{params.organism}" -oe "{params.organelle}" -l "{params.locus}" -b "{params.bioproject}" \
        -s "{params.sample_type}" -t "{params.taxonomy}" -oh "{params.organism_header}" -mt "{params.mol_type}" -top "{params.topology}"
        touch {output}
        """