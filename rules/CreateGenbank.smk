#!/usr/bin/bash
############################# INTRODUCTION #############################
# 
# When the Assembly for the mitochondrial genome and the annotation is generated. Can use the resulting FASTA and GFF file
# for creating the Genbank file. All kind of options are defined in the config which are parsed on the command line to create a Genbank file.
# 
###################################################################
#
############################# MODULES #############################
import os, glob

############################# RULES #############################
# Take out the unclassified reads that will be used to create the Genome Skim database. 
rule CreateGenbank:
    input:
        AnnoMitos = expand(os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
        project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
        ),
        AssemblyGetOrg = expand(os.path.join(DATA_DIR_GS, "{project}/04_GetOrganelle_Mito_Results/{sample}/get_org.log.txt"),
        project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
        )

    output:
        expand(os.path.join(DATA_DIR_GS, "{project}/16_Genbank/{TAXID}_{sample}.done"),
        project = config["genome_skimming"]["project"], 
        sample=config["genome_skimming"]["sample"],
        TAXID = config['genome_skimming']['taxid']
        )

    conda:
        "../envs/biopython.yaml"

    params:
        FastaFile = expand(os.path.join(DATA_DIR_GS, "{project}/04_GetOrganelle_Mito_Results/{sample}/*1.1.path_sequence.fasta"),
        project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
        ),
        GffFile = expand(os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/"),
        project = config["genome_skimming"]["project"], sample=config["genome_skimming"]["sample"]
        ),
        GenbankFile = expand(os.path.join(DATA_DIR_GS, "{project}/16_Genbank/{TAXID}_{sample}.gb"),
        project = config["genome_skimming"]["project"], 
        sample=config["genome_skimming"]["sample"],
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